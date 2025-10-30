import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import numpy as np
import warnings
from scipy import signal

def read_labchart_txt(file_path: str, delimiter: str = '\t') -> pd.DataFrame:
    """
    Reads a LabChart .txt file by manually parsing it line-by-line.
    
    This function is designed to be robust against common LabChart export
    issues, such as metadata headers, mid-file comment blocks, and
    rows with inconsistent column counts.
    
    Parameters:
    - file_path (str): The full path to the .txt file.
    - delimiter (str): The delimiter used in the file (default is tab).
    
    Returns:
    - pd.DataFrame: A cleaned DataFrame containing the time-series data.
    """
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # Find the start of the first block of actual data by skipping metadata
    data_start_index = 0
    for i, line in enumerate(lines):
        stripped_line = line.strip()
        # Assume data starts on the first line that begins with a digit, comma, or minus
        if stripped_line and (stripped_line[0].isdigit() or stripped_line.startswith((',', '-'))):
            data_start_index = i
            break
    if data_start_index == 0:
        warnings.warn(f"Could not detect metadata in {file_path}. Assuming data starts on line 1.")

    # Parse the Header from the initial metadata to get channel names
    channel_names = []
    metadata_lines = lines[:data_start_index]
    header_line = next((l for l in metadata_lines if 'ChannelTitle=' in l), None)
    
    if header_line:
        # Extract channel names from the 'ChannelTitle=' line
        header_text = header_line.split('=', 1)[1].strip()
        channel_names = [name.strip() for name in header_text.split(delimiter) if name.strip()]
    
    # If no header is found, try to guess column count from the first data row
    if not channel_names:
        if data_start_index >= len(lines):
            warnings.warn(f"File {file_path} appears to be empty.")
            return pd.DataFrame() # No data
            
        first_data_row_parts = lines[data_start_index].strip().split(delimiter)
        num_parts = len(first_data_row_parts)
        
        if num_parts <= 1:
             warnings.warn(f"Could not parse header or data in {file_path}.")
             return pd.DataFrame() # Not enough data

        # Guess number of channels from first data row
        if num_parts > 1 and first_data_row_parts[-1].startswith('#*'):
            num_data_cols = num_parts - 2 # Time, Comment
        else:
            num_data_cols = num_parts - 1 # Time, no Comment
            
        if num_data_cols < 0: num_data_cols = 0
        channel_names = [f'Channel_{i+1}' for i in range(num_data_cols)]
        warnings.warn(f"'ChannelTitle=' not found. Guessing {len(channel_names)} channels.")

    # Define the final column structure based on the header
    final_col_names = ['Time'] + channel_names + ['Comment']
    num_expected_cols = len(final_col_names)

    # Process all data rows robustly
    data_rows = []
    for line_num, line in enumerate(lines[data_start_index:], start=data_start_index):
        line = line.strip()
        if not line:
            continue
        
        parts = line.split(delimiter)
        row_len = len(parts)
        
        if row_len == num_expected_cols:
            # Row has the expected number of columns
            data_rows.append(parts)
        
        elif row_len == num_expected_cols - 1:
            # Row is missing the optional 'Comment' column
            parts.append(None)
            data_rows.append(parts)
        
        elif row_len > num_expected_cols:
            # Row has *more* columns than expected.
            # Assume extra parts are a malformed comment and merge them.
            data_part = parts[:num_expected_cols - 1]
            comment_part = " ".join(parts[num_expected_cols - 1:])
            data_rows.append(data_part + [comment_part])
        
        elif row_len < num_expected_cols - 1:
            # Row has *fewer* columns than expected.
            # Pad with None to prevent crashing and allow for later cleanup.
            parts.extend([None] * (num_expected_cols - row_len))
            data_rows.append(parts)

    if not data_rows:
        warnings.warn(f"No data rows found in {file_path}.")
        return pd.DataFrame()

    # Create the DataFrame
    df = pd.DataFrame(data_rows, columns=final_col_names)

    # Convert all columns to numeric, except for 'Comment'
    for col in df.columns:
        if col != 'Comment':
            # Handle LabChart's use of commas for decimal points
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(',', '.', regex=False), errors='coerce')
    df['Comment'] = df['Comment'].astype('string')
    
    # Drop any rows where 'Time' could not be parsed (e.g., junk metadata lines)
    df.dropna(subset=['Time'], inplace=True)
    return df.reset_index(drop=True)


def tms_csp(
    folder_path: str,
    emg_column_name: str,
    output_plot_path: str,
    force_column_name: str = None,
    m_wave_trigger_text: str = None,
    time_column_name: str = 'Time',
    tms_trigger_text: str = '#* TMS Pulse',
    pre_stim_baseline_window_ms: tuple = (-100, 0),
    long_pre_stim_window_s: float = 2.5,
    mep_window_ms: tuple = (5.0, 50.0),
    csp_method: str = 'median',
    csp_sd_multiplier: float = 2.0,
    csp_threshold_percentile: float = 0.50,
    csp_mcd_multiplier: float = 2.66,
    csp_search_window_ms: float = 400.0,
    csp_smoothing_window_ms: float = 10.0,
    csp_min_return_duration_ms: float = 100.0,
    csp_onset_method: str = 'tms_pulse',
    mep_auc_window_ms: tuple = (15.0, 100.0),
    mep_onset_sd_thresh: float = 3.0,
    plot_rectified_emg: bool = False
) -> pd.DataFrame:
    """
    Analyzes all LabChart .txt files in a folder for MEP and CSP parameters.
    
    Processes each file, identifies stimuli (TMS and M-wave), and calculates
    key metrics like MEP amplitude, MEP AUC, and CSP duration. Generates
    a detailed plot for each trial and returns a summary DataFrame.

    Parameters:
    - folder_path (str): Path to the directory containing .txt files.
    - emg_column_name (str): The exact name of the EMG channel to analyze.
    - output_plot_path (str): Path to save the generated plots.
    - force_column_name (str, optional): Name of the force channel.
    - m_wave_trigger_text (str, optional): Comment text to identify M-waves.
    - time_column_name (str): Name of the time column (default 'Time').
    - tms_trigger_text (str): Comment text to identify TMS pulses.
    - pre_stim_baseline_window_ms (tuple): (start, end) ms relative to pulse
      for calculating CSP thresholds and RMS.
    - long_pre_stim_window_s (float): Duration in seconds before pulse
      for calculating long-window RMS and force CV.
    - mep_window_ms (tuple): (start, end) ms post-pulse to find MEP peak-to-peak.
    - csp_method (str): Method for CSP threshold ('median', 'sd', 'mcd').
    - csp_sd_multiplier (float): Multiplier for 'sd' method.
    - csp_threshold_percentile (float): Percentile for 'median' method.
    - csp_mcd_multiplier (float): Multiplier for 'mcd' method (Garvey et al. 2001).
    - csp_search_window_ms (float): Max duration in ms to search for CSP offset.
    - csp_smoothing_window_ms (float): Rolling average window (ms) to smooth
      rectified EMG before checking CSP offset.
    - csp_min_return_duration_ms (float): Duration (ms) the smoothed EMG must
      stay above threshold to confirm CSP offset.
    - csp_onset_method (str): 'tms_pulse', 'mep_onset', or 'mep_offset'.
    - mep_auc_window_ms (tuple): (start, end) ms for MEP AUC calculation.
    - mep_onset_sd_thresh (float): SD multiplier for dynamic MEP onset detection.
    - plot_rectified_emg (bool): If True, plots the rectified EMG under the raw signal.
    
    Returns:
    - pd.DataFrame: A DataFrame with all calculated results, one row per pulse.
    """
    # Ensure the output directory exists
    os.makedirs(output_plot_path, exist_ok=True)
    all_results = []
    file_paths = list(Path(folder_path).glob('*.txt'))
    
    if not file_paths:
        print(f"Warning: No .txt files found in '{folder_path}'")
        return pd.DataFrame()

    print(f"Found {len(file_paths)} files to process...")

    # Iterate over each file in the folder
    for file_path in file_paths:
        filename_no_ext = file_path.stem
        print(f"\nProcessing: {file_path.name}")

        try:
            # Read the file using the robust reader
            df = read_labchart_txt(file_path)
            local_force_column = force_column_name
            
            # Basic validation
            if df.empty or time_column_name not in df.columns or emg_column_name not in df.columns:
                 print(f"  - Error: Required columns ('{time_column_name}', '{emg_column_name}') not found. Skipping file.")
                 continue
            if local_force_column and local_force_column not in df.columns:
                print(f"  - Warning: Force column '{local_force_column}' not found. Proceeding without force analysis.")
                local_force_column = None
        except Exception as e:
            print(f"  - Error reading {file_path.name}: {e}. Skipping file.")
            continue
        
        # --- Global M-Wave Search (Once per file) ---
        # Search the *entire* file for M-waves to find the true M-max
        all_m_waves_in_file = []
        if m_wave_trigger_text:
            print("  - Searching for all M-waves in the file...")
            m_wave_indices = df[df['Comment'].astype(str).str.contains(m_wave_trigger_text, na=False)].index
            m_wave_times = df.loc[m_wave_indices, time_column_name].tolist()
            
            if m_wave_times:
                for pulse_time in m_wave_times:
                    # Define M-wave window (0 to 50ms post-pulse)
                    m_wave_df = df[(df[time_column_name] >= pulse_time) & (df[time_column_name] <= pulse_time + 0.050)]
                    amplitude = 0
                    m_wave_auc = 0
                    
                    if not m_wave_df.empty:
                        # M-wave amplitude (peak-to-peak)
                        amplitude = m_wave_df[emg_column_name].max() - m_wave_df[emg_column_name].min()
                        
                        # M-wave AUC (baseline-corrected)
                        m_wave_baseline_df = df[(df[time_column_name] >= pulse_time - 0.100) & (df[time_column_name] < pulse_time)]
                        if not m_wave_baseline_df.empty:
                            baseline_mean_rect = m_wave_baseline_df[emg_column_name].abs().mean()
                            rectified_m_wave = m_wave_df[emg_column_name].abs()
                            corrected_rectified = rectified_m_wave - baseline_mean_rect
                            corrected_rectified[corrected_rectified < 0] = 0
                            m_wave_auc = np.trapz(y=corrected_rectified, x=m_wave_df[time_column_name])

                    all_m_waves_in_file.append({'time': pulse_time, 'amplitude': amplitude, 'auc': m_wave_auc})
                print(f"  - Found {len(all_m_waves_in_file)} M-wave(s) in total for this file.")
            else:
                print("  - No M-wave stimuli found anywhere in this file.")

        # --- Section Splitting ---
        # Split the file into continuous blocks based on time discontinuities
        first_valid_index = df.index[0]
        time_diffs = df[time_column_name].diff()
        median_step = time_diffs.median()
        
        if pd.isna(median_step) or median_step <= 0:
            # If time step is invalid, treat as one big section
            section_start_indices = [first_valid_index]
        else:
            # A discontinuity is a time jump > 100x the median step
            discontinuity_threshold = median_step * 100
            discontinuities = list(time_diffs[(time_diffs < 0) | (time_diffs > discontinuity_threshold)].index)
            section_start_indices = [first_valid_index] + discontinuities
        
        section_start_indices = sorted(list(set(section_start_indices)))
        sections = [df.loc[start:end-1] if end is not None else df.loc[start:] for start, end in zip(section_start_indices, section_start_indices[1:] + [None])]
        
        print(f"  - Detected {len(sections)} section(s) in this file.")
        
        # --- Main Analysis Loop (Per Section) ---
        for section_num, section_df in enumerate(sections, 1):
            if section_df.empty:
                continue
            
            section_df = section_df.copy()
            # Pre-calculate rectified EMG for this section
            rectified_emg_col = f"{emg_column_name}_rectified"
            section_df[rectified_emg_col] = section_df[emg_column_name].abs()

            # Calculate sampling rate for this section
            time_diffs_section = section_df[time_column_name].diff()
            median_sample_period_s = time_diffs_section.median()
            
            smoothing_samples_needed = 1
            duration_samples_needed = 1
            
            if pd.isna(median_sample_period_s) or median_sample_period_s <= 0:
                print(f"  - Warning: Could not determine valid sampling rate for Section #{section_num}. CSP calcs will use 1 sample.")
            else:
                # Convert time-based windows (ms) to sample-based windows
                smoothing_duration_s = csp_smoothing_window_ms / 1000.0
                smoothing_samples_needed = int(round(smoothing_duration_s / median_sample_period_s))
                if smoothing_samples_needed < 1:
                    smoothing_samples_needed = 1
                    
                duration_s = csp_min_return_duration_ms / 1000.0
                duration_samples_needed = int(round(duration_s / median_sample_period_s))
                if duration_samples_needed < 1:
                    duration_samples_needed = 1

            # Find all M-Wave and TMS stimuli in this section
            all_stimuli = []
            if m_wave_trigger_text:
                m_wave_indices = section_df[section_df['Comment'].astype(str).str.contains(m_wave_trigger_text, na=False)].index
                for index in m_wave_indices:
                    all_stimuli.append({'time': section_df.loc[index, time_column_name], 'type': 'M-Wave'})
            
            tms_indices = section_df[section_df['Comment'].astype(str).str.contains(tms_trigger_text, na=False)].index
            for index in tms_indices:
                all_stimuli.append({'time': section_df.loc[index, time_column_name], 'type': 'TMS'})

            if not all_stimuli:
                print(f"  - No stimuli found in Section #{section_num}.")
                continue
            
            # Sort stimuli by time to process in order
            all_stimuli.sort(key=lambda x: x['time'])
            
            # --- Intelligent Mmax determination for this section ---
            m_max_amp_for_section = None
            m_max_auc_for_section = None
            if all_m_waves_in_file:
                section_start_time = section_df[time_column_name].iloc[0]
                section_end_time = section_df[time_column_name].iloc[-1]
                
                # Try to find M-waves *within* this section first
                m_waves_in_section = [mw for mw in all_m_waves_in_file if section_start_time <= mw['time'] <= section_end_time]
                if m_waves_in_section:
                    m_max_amp_for_section = max(mw['amplitude'] for mw in m_waves_in_section)
                    m_max_auc_for_section = max(mw['auc'] for mw in m_waves_in_section)
                    print(f"  - Section #{section_num}: Using section-specific Mmax Amp {m_max_amp_for_section:.4f} and Mmax AUC {m_max_auc_for_section:.4f}")
                else:
                    # If no M-wave in this section, use the closest one from the file
                    section_mean_time = section_df[time_column_name].mean()
                    closest_m_wave = min(all_m_waves_in_file, key=lambda mw: abs(mw['time'] - section_mean_time))
                    m_max_amp_for_section = closest_m_wave['amplitude']
                    m_max_auc_for_section = closest_m_wave['auc']
                    print(f"  - Section #{section_num}: No M-wave in section. Using closest M-wave's Amp {m_max_amp_for_section:.4f} and AUC {m_max_auc_for_section:.4f}")
            else:
                 print(f"  - Section #{section_num}: No M-waves in file to use for normalization.")

            print(f"  - Section #{section_num}: Found {len(all_stimuli)} total stimuli. Analyzing and plotting...")
            
            # --- Plot Initialization (One figure per section) ---
            plots_per_pulse = 2 if local_force_column else 1
            num_rows = len(all_stimuli) * plots_per_pulse
            fig_height = (3 * len(all_stimuli)) * plots_per_pulse
            fig, axes = plt.subplots(nrows=num_rows, ncols=1, figsize=(12, fig_height), squeeze=False)
            axes = axes.flatten()

            # --- Per-Stimulus Analysis and Plotting Loop ---
            for i, stim in enumerate(all_stimuli):
                pulse_time = stim['time']
                stim_type = stim['type']
                
                ax_idx = i * plots_per_pulse
                ax_emg = axes[ax_idx]
                # Define a standard plotting window around the pulse
                plot_window_df = section_df[(section_df[time_column_name] >= pulse_time - 0.2) & (section_df[time_column_name] <= pulse_time + 0.5)]

                # --- M-Wave Processing ---
                if stim_type == 'M-Wave':
                    m_wave_df = section_df[(section_df[time_column_name] >= pulse_time) & (section_df[time_column_name] <= pulse_time + 0.050)]
                    amplitude = m_wave_df[emg_column_name].max() - m_wave_df[emg_column_name].min() if not m_wave_df.empty else 0
                    
                    # Calculate M-wave AUC (same logic as global search)
                    m_wave_auc = 0
                    m_wave_baseline_df = section_df[(section_df[time_column_name] >= pulse_time - 0.100) & (section_df[time_column_name] < pulse_time)]
                    if not m_wave_baseline_df.empty:
                        baseline_mean_rect = m_wave_baseline_df[emg_column_name].abs().mean()
                        rectified_m_wave = m_wave_df[emg_column_name].abs()
                        corrected_rectified = rectified_m_wave - baseline_mean_rect
                        corrected_rectified[corrected_rectified < 0] = 0
                        m_wave_auc = np.trapz(y=corrected_rectified, x=m_wave_df[time_column_name])
                    
                    # Plot M-Wave Data
                    if plot_rectified_emg:
                        ax_emg.plot(plot_window_df[time_column_name], plot_window_df[rectified_emg_col], color='lightgrey', lw=0.5, label='_nolegend_')
                    
                    ax_emg.plot(plot_window_df[time_column_name], plot_window_df[emg_column_name], color='black', lw=0.5)
                    ax_emg.axvline(pulse_time, color='blue', linestyle='--', label='M-Wave Stimulus')
                    ax_emg.fill_between(m_wave_df[time_column_name], 0, m_wave_df[emg_column_name], color='purple', alpha=0.3, label=f'M-Wave AUC: {m_wave_auc:.4f}')
                    if not m_wave_df.empty:
                        max_idx, min_idx = m_wave_df[emg_column_name].idxmax(), m_wave_df[emg_column_name].idxmin()
                        ax_emg.plot(m_wave_df.loc[max_idx, time_column_name], m_wave_df.loc[max_idx, emg_column_name], 'o', color='blue', markersize=5)
                        ax_emg.plot(m_wave_df.loc[min_idx, time_column_name], m_wave_df.loc[min_idx, emg_column_name], 'o', color='red', markersize=5, label=f'M-Wave Amp: {amplitude:.2f}')
                    ax_emg.set_title(f"M-Wave at {pulse_time:.3f} s", fontsize=10)

                # --- TMS Pulse Processing ---
                elif stim_type == 'TMS':
                    # 1. MEP Amplitude Calculation
                    mep_start_time = pulse_time + (mep_window_ms[0] / 1000.0)
                    mep_end_time = pulse_time + (mep_window_ms[1] / 1000.0)
                    mep_df = section_df[(section_df[time_column_name] >= mep_start_time) & (section_df[time_column_name] <= mep_end_time)]
                    mep_amplitude = 0
                    mep_min_val, mep_max_val, mep_min_time, mep_max_time = None, None, None, None
                    if not mep_df.empty:
                        mep_amplitude = mep_df[emg_column_name].max() - mep_df[emg_column_name].min()
                        max_idx, min_idx = mep_df[emg_column_name].idxmax(), mep_df[emg_column_name].idxmin()
                        mep_max_val, mep_min_val = mep_df.loc[max_idx, emg_column_name], mep_df.loc[min_idx, emg_column_name]
                        mep_max_time, mep_min_time = mep_df.loc[max_idx, time_column_name], mep_df.loc[min_idx, time_column_name]

                    # 2. Short Pre-stimulus Baseline (for thresholds and RMS)
                    short_pre_stim_df = section_df[(section_df[time_column_name] >= pulse_time + (pre_stim_baseline_window_ms[0] / 1000.0)) & (section_df[time_column_name] < pulse_time + (pre_stim_baseline_window_ms[1] / 1000.0))]
                    short_pre_stim_emg_rms = float('nan')
                    short_pre_stim_mean_force = float('nan')
                    if short_pre_stim_df.empty:
                        print(f"    - Warning: No short pre-stimulus data for pulse at {pulse_time:.3f}s.")
                    else:
                        short_pre_stim_emg_rms = np.sqrt(np.mean(short_pre_stim_df[emg_column_name]**2))
                        if local_force_column:
                            short_pre_stim_mean_force = short_pre_stim_df[local_force_column].mean()

                    # 3. Long Pre-stimulus Baseline (for stability measures)
                    long_pre_stim_df = section_df[(section_df[time_column_name] >= pulse_time - long_pre_stim_window_s) & (section_df[time_column_name] < pulse_time)]
                    long_pre_stim_emg_rms = float('nan')
                    long_pre_stim_mean_force = float('nan')
                    long_pre_stim_force_cv = float('nan')
                    if long_pre_stim_df.empty:
                        print(f"    - Warning: No long pre-stimulus data for pulse at {pulse_time:.3f}s.")
                    else:
                        long_pre_stim_emg_rms = np.sqrt(np.mean(long_pre_stim_df[emg_column_name]**2))
                        if local_force_column:
                            force_data = long_pre_stim_df[local_force_column]
                            long_pre_stim_mean_force = force_data.mean()
                            if long_pre_stim_mean_force != 0:
                               long_pre_stim_force_cv = (force_data.std() / abs(long_pre_stim_mean_force)) * 100
                    
                    # 4. Dynamic MEP Window and AUC Calculation
                    mep_auc = 0
                    auc_start_time, auc_end_time = None, None # Initialize
                    if not short_pre_stim_df.empty:
                        rectified_baseline_emg = short_pre_stim_df[emg_column_name].abs()
                        pre_stim_mean_rectified = rectified_baseline_emg.mean()
                        
                        # Find dynamic onset/offset based on SD threshold
                        baseline_mean = short_pre_stim_df[emg_column_name].mean()
                        baseline_sd = short_pre_stim_df[emg_column_name].std()
                        upper_thresh = baseline_mean + mep_onset_sd_thresh * baseline_sd
                        lower_thresh = baseline_mean - mep_onset_sd_thresh * baseline_sd
                        
                        onset_search_df = section_df[(section_df[time_column_name] >= pulse_time + 0.020) & (section_df[time_column_name] <= pulse_time + 0.050)]
                        onset_point = onset_search_df[(onset_search_df[emg_column_name] > upper_thresh) | (onset_search_df[emg_column_name] < lower_thresh)]
                        auc_start_time = onset_point.iloc[0][time_column_name] if not onset_point.empty else None
                        
                        offset_search_df = section_df[(section_df[time_column_name] >= pulse_time + 0.050) & (section_df[time_column_name] <= pulse_time + 0.150)]
                        offset_point = offset_search_df[(offset_search_df[emg_column_name] < upper_thresh) & (offset_search_df[emg_column_name] > lower_thresh)]
                        auc_end_time = offset_point.iloc[0][time_column_name] if not offset_point.empty else None
                        
                        # Fallback to fixed window if dynamic detection fails
                        if auc_start_time is None or auc_end_time is None:
                            auc_start_time = pulse_time + (mep_auc_window_ms[0] / 1000.0)
                            auc_end_time = pulse_time + (mep_auc_window_ms[1] / 1000.0)
                        
                        mep_auc_df = section_df[(section_df[time_column_name] >= auc_start_time) & (section_df[time_column_name] <= auc_end_time)].copy()
                        if not mep_auc_df.empty:
                            rectified_mep = mep_auc_df[emg_column_name].abs()
                            corrected_rectified_mep = rectified_mep - pre_stim_mean_rectified
                            corrected_rectified_mep[corrected_rectified_mep < 0] = 0
                            mep_auc = np.trapz(y=corrected_rectified_mep, x=mep_auc_df[time_column_name])
                    else:
                        mep_auc_df = pd.DataFrame()
                        # Ensure fallback if baseline is empty
                        auc_start_time = pulse_time + (mep_auc_window_ms[0] / 1000.0)
                        auc_end_time = pulse_time + (mep_auc_window_ms[1] / 1000.0)

                    # 5. Define CSP Onset Time based on user method
                    csp_start_time = pulse_time # Default to 'tms_pulse'
                    if csp_onset_method.lower() == 'mep_onset':
                        csp_start_time = auc_start_time
                    elif csp_onset_method.lower() == 'mep_offset':
                        csp_start_time = auc_end_time
                    elif csp_onset_method.lower() != 'tms_pulse':
                        warnings.warn(f"Invalid csp_onset_method '{csp_onset_method}'. Defaulting to 'tms_pulse'.")

                    # 6. CSP Detection
                    csp_end_time, csp_duration_ms, csp_threshold_value = None, float('nan'), 0
                    if not short_pre_stim_df.empty:
                        rectified_baseline_emg = short_pre_stim_df[emg_column_name].abs()

                        # 6a. Calculate CSP Threshold based on user method
                        if csp_method.lower() == 'sd':
                            baseline_mean = rectified_baseline_emg.mean()
                            baseline_sd = rectified_baseline_emg.std()
                            csp_threshold_value = baseline_mean + (csp_sd_multiplier * baseline_sd)
                        
                        elif csp_method.lower() == 'median':
                            csp_threshold_value = rectified_baseline_emg.quantile(csp_threshold_percentile)
                        
                        elif csp_method.lower() == 'mcd':
                            baseline_mean = rectified_baseline_emg.mean()
                            consecutive_diffs = rectified_baseline_emg.diff()
                            abs_consecutive_diffs = consecutive_diffs.abs()
                            mcd = abs_consecutive_diffs.mean()
                            # Use baseline_mean - ... for the *lower* limit, per Garvey
                            csp_threshold_value = baseline_mean - (mcd * csp_mcd_multiplier) 
                            if csp_threshold_value < 0:
                                csp_threshold_value = 0 # Guardrail for negative threshold
                        
                        else:
                            warnings.warn(f"Invalid csp_method '{csp_method}'. Defaulting to 'median'.")
                            csp_threshold_value = rectified_baseline_emg.quantile(csp_threshold_percentile)
                        
                        
                        # 6b. Find CSP Offset (Double-threshold logic)
                        # Define search window (e.g., 50ms post-pulse to 400ms post-pulse)
                        csp_search_start = pulse_time + (mep_window_ms[1] / 1000.0) + 0.05 # Start search 50ms after pulse
                        csp_search_end = pulse_time + (csp_search_window_ms / 1000.0)
                        csp_search_df = section_df[(section_df[time_column_name] > csp_search_start) & (section_df[time_column_name] <= csp_search_end)].copy()
                        
                        if not csp_search_df.empty:
                            # Step 1: Smooth the rectified signal
                            smoothed_emg = csp_search_df[rectified_emg_col].rolling(window=smoothing_samples_needed).mean()
                            
                            # Step 2: Find where smoothed signal is above threshold
                            above_threshold = (smoothed_emg > csp_threshold_value)
                            
                            # Step 3: Find where it stays above for the required *duration*
                            consecutive_above = above_threshold.rolling(window=duration_samples_needed).sum()
                            
                            # Find the *first index* where this condition is met
                            first_window_end_index = consecutive_above[consecutive_above >= duration_samples_needed].first_valid_index()

                            if first_window_end_index is not None:
                                # Found the end of the first valid window.
                                # Now find the *start* of this window.
                                end_position = csp_search_df.index.get_loc(first_window_end_index)
                                start_position = end_position - (duration_samples_needed - 1)
                                start_index = csp_search_df.index[start_position]
                                
                                # This is the CSP offset time
                                csp_end_time = csp_search_df.loc[start_index, time_column_name]
                                
                                # Calculate final duration from the chosen start time
                                csp_duration_ms = (csp_end_time - csp_start_time) * 1000.0
                        
                    # 7. MEP Normalization
                    mep_amplitude_normalized, mep_auc_normalized = float('nan'), float('nan')
                    if m_max_amp_for_section and m_max_amp_for_section > 0:
                        mep_amplitude_normalized = (mep_amplitude / m_max_amp_for_section) * 100
                    if m_max_auc_for_section and m_max_auc_for_section > 0:
                        mep_auc_normalized = (mep_auc / m_max_auc_for_section) * 100

                    # 8. Store Results
                    results_to_append = {
                        'filename': file_path.name, 'section': section_num, 'pulse_time_s': pulse_time,
                        'mep_amplitude': mep_amplitude, 'mep_auc': mep_auc, 'csp_duration_ms': csp_duration_ms,
                        'mep_amplitude_percent_mmax': mep_amplitude_normalized, 'mep_auc_percent_mmax': mep_auc_normalized,
                        'mmax_amp_used': m_max_amp_for_section, 'mmax_auc_used': m_max_auc_for_section,
                        'short_pre_stim_emg_rms': short_pre_stim_emg_rms,
                        'long_pre_stim_emg_rms': long_pre_stim_emg_rms,
                    }
                    if local_force_column:
                        results_to_append.update({
                            'short_pre_stim_mean_force': short_pre_stim_mean_force,
                            'long_pre_stim_mean_force': long_pre_stim_mean_force,
                            'long_pre_stim_force_cv': long_pre_stim_force_cv
                        })
                    
                    all_results.append(results_to_append)

                    # 9. Facetted Plotting Logic
                    
                    # Plot grey CSP box *first* so it's in the background
                    if csp_end_time:
                        ax_emg.axvspan(csp_start_time, csp_end_time, color='gray', alpha=0.3, label=f'CSP: {csp_duration_ms:.1f} ms')
                        ax_emg.axvline(csp_end_time, color='orange', linestyle='--', label=f'CSP End')

                    # Plot optional rectified signal
                    if plot_rectified_emg:
                        ax_emg.plot(plot_window_df[time_column_name], plot_window_df[rectified_emg_col], color='lightgrey', lw=0.5, label='_nolegend_')
                    
                    # Plot raw EMG signal
                    ax_emg.plot(plot_window_df[time_column_name], plot_window_df[emg_column_name], color='black', lw=0.5)
                    
                    # Plot vertical/horizontal lines and markers
                    ax_emg.axvline(pulse_time, color='red', linestyle='--', label=f'TMS Pulse')
                    ax_emg.axhline(csp_threshold_value, color='blue', linestyle=':', lw=1, label=f'CSP Threshold')
                    
                    # Plot MEP AUC fill
                    if not mep_auc_df.empty:
                        ax_emg.fill_between(mep_auc_df[time_column_name], 0, mep_auc_df[emg_column_name], color='green', alpha=0.3, label=f'MEP AUC: {mep_auc:.4f}')
                    
                    # Plot MEP amplitude markers
                    if mep_max_time is not None:
                        ax_emg.plot(mep_max_time, mep_max_val, 'o', color='blue', markersize=5)
                        ax_emg.plot(mep_min_time, mep_min_val, 'o', color='red', markersize=5, label=f'MEP Amp: {mep_amplitude:.2f}')
                    
                    ax_emg.set_title(f"TMS Pulse at {pulse_time:.3f} s", fontsize=10)

                # --- Common plotting setup for both M-Wave and TMS ---
                ax_emg.set_ylabel(f"EMG Amp ({emg_column_name})")
                ax_emg.legend(loc='upper right', fontsize=8)
                ax_emg.grid(True, linestyle='--', alpha=0.6)

                # Plot force channel if it exists
                if local_force_column:
                    ax_force = axes[ax_idx + 1]
                    ax_emg.get_shared_x_axes().joined(ax_emg, ax_force)
                    ax_force.plot(plot_window_df[time_column_name], plot_window_df[local_force_column], color='green', lw=0.5, linestyle='--')
                    ax_force.set_ylabel(f"Force ({local_force_column})")
                    ax_force.set_xlabel("Time (s)")
                    ax_force.grid(True, linestyle='--', alpha=0.6)
                    # Hide x-axis labels on the top EMG plot
                    ax_emg.tick_params(axis='x', labelbottom=False)
                else:
                    ax_emg.set_xlabel("Time (s)")
            
            # --- Save Figure ---
            fig.suptitle(f'Analysis for: {file_path.name} - Section #{section_num}', fontsize=16, y=1.0)
            fig.tight_layout(rect=[0, 0, 1, 0.98], h_pad=3.0)
            plot_filename = f"{filename_no_ext}_Section_{section_num}.jpg"
            fig.savefig(os.path.join(output_plot_path, plot_filename), dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"  - Combined plot for Section #{section_num} saved to {plot_filename}")

    print("\n✅ Analysis complete.")
    return pd.DataFrame(all_results)
