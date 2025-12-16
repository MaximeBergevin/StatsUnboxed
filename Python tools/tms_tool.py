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
    """
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # Find the start of the first block of actual data
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
    csp_sd_multiplier: float = 2.0,
    csp_mcd_multiplier: float = 2.66,
    csp_search_window_ms: float = 400.0,
    csp_smoothing_window_ms: float = 10.0,
    csp_min_return_duration_ms: float = 100.0,
    mep_auc_window_ms: tuple = (15.0, 100.0),
    mep_onset_sd_thresh: float = 3.0,
    plot_rectified_emg: bool = False
) -> pd.DataFrame:
    """
    Analyzes LabChart .txt files with resetting timestamps (Scope mode), 
    ensuring M-waves are paired with the correct section's TMS pulses.
    """
    os.makedirs(output_plot_path, exist_ok=True)
    all_results = []
    file_paths = list(Path(folder_path).glob('*.txt'))
    
    if not file_paths:
        print(f"Warning: No .txt files found in '{folder_path}'")
        return pd.DataFrame()

    print(f"Found {len(file_paths)} files to process...")

    for file_path in file_paths:
        filename_no_ext = file_path.stem
        print(f"\nProcessing: {file_path.name}")

        try:
            # Assumes read_labchart_txt is defined elsewhere in your script
            df = read_labchart_txt(file_path)
            
            local_force_column = force_column_name
            if df.empty or time_column_name not in df.columns or emg_column_name not in df.columns:
                 print(f"  - Error: Required columns ('{time_column_name}', '{emg_column_name}') not found. Skipping file.")
                 continue
            if local_force_column and local_force_column not in df.columns:
                print(f"  - Warning: Force column '{local_force_column}' not found. Proceeding without force analysis.")
                local_force_column = None
        except Exception as e:
            print(f"  - Error reading {file_path.name}: {e}. Skipping file.")
            continue

        # =================================================================
        # STEP 1: Split Sections Based on Time Discontinuities
        # =================================================================
        # This is critical for resetting timestamps.
        first_valid_index = df.index[0]
        time_diffs = df[time_column_name].diff()
        median_step = time_diffs.median()
        
        if pd.isna(median_step) or median_step <= 0:
            section_start_indices = [first_valid_index]
        else:
            # If time jumps backwards or jumps forward massively, it's a new section
            discontinuity_threshold = median_step * 100
            discontinuities = list(time_diffs[(time_diffs < 0) | (time_diffs > discontinuity_threshold)].index)
            section_start_indices = [first_valid_index] + discontinuities
        
        section_start_indices = sorted(list(set(section_start_indices)))
        
        # Create list of dataframes, one per section
        sections = []
        for start, end in zip(section_start_indices, section_start_indices[1:] + [None]):
            if end is None:
                sections.append(df.loc[start:])
            else:
                sections.append(df.loc[start:end-1])
        
        print(f"  - Detected {len(sections)} section(s) in this file.")

        # =================================================================
        # STEP 2: Pre-calculate M-waves PER SECTION
        # =================================================================
        # Dictionary to store M-max data: { section_index (int): { 'amplitude': float, 'auc': float } }
        section_m_waves = {}
        valid_m_wave_sections = [] # To help with fallback logic

        if m_wave_trigger_text:
            print("  - Scanning sections for local M-waves...")
            for s_idx, s_df in enumerate(sections):
                # Find triggers ONLY in this section
                m_indices = s_df[s_df['Comment'].astype(str).str.contains(m_wave_trigger_text, na=False)].index
                
                # If triggers exist, process the LAST one in the section (common protocol standard)
                if not m_indices.empty:
                    pulse_time = s_df.loc[m_indices[-1], time_column_name]
                    
                    # Define M-wave window (0 to +50ms)
                    mw_df = s_df[(s_df[time_column_name] >= pulse_time) & (s_df[time_column_name] <= pulse_time + 0.050)]
                    
                    if not mw_df.empty:
                        # 1. Amplitude
                        amp = mw_df[emg_column_name].max() - mw_df[emg_column_name].min()
                        
                        # 2. AUC (Corrected)
                        auc = 0
                        base_df = s_df[(s_df[time_column_name] >= pulse_time - 0.100) & (s_df[time_column_name] < pulse_time)]
                        if not base_df.empty:
                            base_mean = base_df[emg_column_name].abs().mean()
                            rect_sig = mw_df[emg_column_name].abs()
                            corr_sig = rect_sig - base_mean
                            corr_sig[corr_sig < 0] = 0 # Rectified cannot be negative
                            auc = np.trapz(y=corr_sig, x=mw_df[time_column_name])
                        
                        section_m_waves[s_idx] = {'amplitude': amp, 'auc': auc, 'time': pulse_time}
                        valid_m_wave_sections.append(s_idx)
            
            print(f"  - Found M-waves in {len(valid_m_wave_sections)} of {len(sections)} sections.")

        # =================================================================
        # STEP 3: Main Analysis Loop (Section by Section)
        # =================================================================
        for section_num, section_df in enumerate(sections, 1):
            s_idx = section_num - 1
            if section_df.empty:
                continue
            
            section_df = section_df.copy()
            rectified_emg_col = f"{emg_column_name}_rectified"
            section_df[rectified_emg_col] = section_df[emg_column_name].abs()

            # Sampling rate calculation for this section
            time_diffs_section = section_df[time_column_name].diff()
            median_sample_period_s = time_diffs_section.median()
            
            smoothing_samples_needed = 1
            duration_samples_needed = 1
            
            if pd.isna(median_sample_period_s) or median_sample_period_s <= 0:
                print(f"  - Warning: Bad sampling rate in Section #{section_num}. Using default 1 sample.")
            else:
                smoothing_duration_s = csp_smoothing_window_ms / 1000.0
                smoothing_samples_needed = int(round(smoothing_duration_s / median_sample_period_s)) or 1
                duration_s = csp_min_return_duration_ms / 1000.0
                duration_samples_needed = int(round(duration_s / median_sample_period_s)) or 1

            # Find Stimuli in this section
            all_stimuli = []
            if m_wave_trigger_text:
                # Add M-waves to plot list
                m_indices = section_df[section_df['Comment'].astype(str).str.contains(m_wave_trigger_text, na=False)].index
                for index in m_indices:
                    all_stimuli.append({'time': section_df.loc[index, time_column_name], 'type': 'M-Wave'})
            
            tms_indices = section_df[section_df['Comment'].astype(str).str.contains(tms_trigger_text, na=False)].index
            for index in tms_indices:
                all_stimuli.append({'time': section_df.loc[index, time_column_name], 'type': 'TMS'})

            if not all_stimuli:
                print(f"  - No stimuli found in Section #{section_num}.")
                continue
            
            all_stimuli.sort(key=lambda x: x['time'])
            
            # --- DETERMINE M-MAX FOR THIS SECTION ---
            m_max_amp_used = None
            m_max_auc_used = None
            
            # Priority 1: Local M-wave
            if s_idx in section_m_waves:
                m_max_amp_used = section_m_waves[s_idx]['amplitude']
                m_max_auc_used = section_m_waves[s_idx]['auc']
                print(f"  - Section #{section_num}: Using LOCAL M-max Amp {m_max_amp_used:.4f}")
            # Priority 2: Nearest Section Fallback
            elif valid_m_wave_sections:
                closest_s_idx = min(valid_m_wave_sections, key=lambda x: abs(x - s_idx))
                m_max_amp_used = section_m_waves[closest_s_idx]['amplitude']
                m_max_auc_used = section_m_waves[closest_s_idx]['auc']
                print(f"  - Section #{section_num}: No local M-wave. Using Section #{closest_s_idx+1} (Amp {m_max_amp_used:.4f})")
            else:
                print(f"  - Section #{section_num}: No M-waves found in file.")

            # --- Plot Initialization ---
            plots_per_pulse = 2 if local_force_column else 1
            num_rows = len(all_stimuli) * plots_per_pulse
            fig_height = max(4, (3 * len(all_stimuli)) * plots_per_pulse) # Ensure min height
            fig, axes = plt.subplots(nrows=num_rows, ncols=1, figsize=(12, fig_height), squeeze=False)
            axes = axes.flatten()

            # --- Pulse Analysis Loop ---
            for i, stim in enumerate(all_stimuli):
                pulse_time = stim['time']
                stim_type = stim['type']
                ax_idx = i * plots_per_pulse
                ax_emg = axes[ax_idx]
                
                plot_window_df = section_df[(section_df[time_column_name] >= pulse_time - 0.2) & (section_df[time_column_name] <= pulse_time + 0.5)]

                # -------------------------
                # TYPE: M-WAVE (Just Plot)
                # -------------------------
                if stim_type == 'M-Wave':
                    ax_emg.plot(plot_window_df[time_column_name], plot_window_df[emg_column_name], color='black', lw=0.5)
                    ax_emg.axvline(pulse_time, color='blue', linestyle='--', label='M-Wave Stim')
                    
                    # Highlight analysis window
                    mw_local = section_df[(section_df[time_column_name] >= pulse_time) & (section_df[time_column_name] <= pulse_time + 0.050)]
                    if not mw_local.empty:
                        ax_emg.fill_between(mw_local[time_column_name], 0, mw_local[emg_column_name], color='purple', alpha=0.3)
                        amp_local = mw_local[emg_column_name].max() - mw_local[emg_column_name].min()
                        ax_emg.set_title(f"M-Wave at {pulse_time:.3f}s (Amp: {amp_local:.2f})", fontsize=10)
                    else:
                        ax_emg.set_title(f"M-Wave at {pulse_time:.3f}s", fontsize=10)

                # -------------------------
                # TYPE: TMS (Analyze)
                # -------------------------
                elif stim_type == 'TMS':
                    # A. MEP Amplitude
                    mep_start = pulse_time + (mep_window_ms[0] / 1000.0)
                    mep_end = pulse_time + (mep_window_ms[1] / 1000.0)
                    mep_df = section_df[(section_df[time_column_name] >= mep_start) & (section_df[time_column_name] <= mep_end)]
                    
                    mep_amplitude = 0
                    mep_max_time, mep_min_time = None, None
                    mep_max_val, mep_min_val = None, None
                    
                    if not mep_df.empty:
                        mep_amplitude = mep_df[emg_column_name].max() - mep_df[emg_column_name].min()
                        max_idx, min_idx = mep_df[emg_column_name].idxmax(), mep_df[emg_column_name].idxmin()
                        mep_max_val, mep_min_val = mep_df.loc[max_idx, emg_column_name], mep_df.loc[min_idx, emg_column_name]
                        mep_max_time, mep_min_time = mep_df.loc[max_idx, time_column_name], mep_df.loc[min_idx, time_column_name]

                    # B. Baselines (Short & Long)
                    short_pre_start = pulse_time + (pre_stim_baseline_window_ms[0]/1000.0)
                    short_pre_end = pulse_time + (pre_stim_baseline_window_ms[1]/1000.0)
                    short_pre_df = section_df[(section_df[time_column_name] >= short_pre_start) & (section_df[time_column_name] < short_pre_end)]
                    
                    short_rms = np.nan
                    short_mean_force = np.nan
                    if not short_pre_df.empty:
                        short_rms = np.sqrt(np.mean(short_pre_df[emg_column_name]**2))
                        if local_force_column:
                            short_mean_force = short_pre_df[local_force_column].mean()
                            
                    long_pre_df = section_df[(section_df[time_column_name] >= pulse_time - long_pre_stim_window_s) & (section_df[time_column_name] < pulse_time)]
                    long_rms = np.nan
                    long_mean_force, long_force_cv = np.nan, np.nan
                    if not long_pre_df.empty:
                        long_rms = np.sqrt(np.mean(long_pre_df[emg_column_name]**2))
                        if local_force_column:
                            f_data = long_pre_df[local_force_column]
                            long_mean_force = f_data.mean()
                            if abs(long_mean_force) > 1e-6:
                                long_force_cv = (f_data.std() / abs(long_mean_force)) * 100

                    # C. MEP AUC (Dynamic Window)
                    mep_auc = 0
                    auc_start, auc_end = None, None
                    if not short_pre_df.empty:
                        base_mean = short_pre_df[emg_column_name].mean()
                        base_sd = short_pre_df[emg_column_name].std()
                        upper = base_mean + mep_onset_sd_thresh * base_sd
                        lower = base_mean - mep_onset_sd_thresh * base_sd
                        
                        # Find Onset (20-50ms)
                        onset_search = section_df[(section_df[time_column_name] >= pulse_time + 0.020) & (section_df[time_column_name] <= pulse_time + 0.050)]
                        onset_pts = onset_search[(onset_search[emg_column_name] > upper) | (onset_search[emg_column_name] < lower)]
                        auc_start = onset_pts.iloc[0][time_column_name] if not onset_pts.empty else None
                        
                        # Find Offset (50-150ms)
                        if auc_start:
                             offset_search = section_df[(section_df[time_column_name] >= pulse_time + 0.050) & (section_df[time_column_name] <= pulse_time + 0.150)]
                             # Offset is return to baseline
                             offset_pts = offset_search[(offset_search[emg_column_name] < upper) & (offset_search[emg_column_name] > lower)]
                             auc_end = offset_pts.iloc[0][time_column_name] if not offset_pts.empty else None

                        # Fallback
                        if auc_start is None or auc_end is None:
                            auc_start = pulse_time + (mep_auc_window_ms[0]/1000.0)
                            auc_end = pulse_time + (mep_auc_window_ms[1]/1000.0)
                            
                        # Calculate
                        auc_df = section_df[(section_df[time_column_name] >= auc_start) & (section_df[time_column_name] <= auc_end)]
                        if not auc_df.empty:
                             rect_base = short_pre_df[emg_column_name].abs().mean()
                             rect_mep = auc_df[emg_column_name].abs()
                             corr_mep = rect_mep - rect_base
                             corr_mep[corr_mep < 0] = 0
                             mep_auc = np.trapz(y=corr_mep, x=auc_df[time_column_name])

                    # D. CSP Calculations
                    # Onset definitions
                    csp_onset_pulse = pulse_time
                    csp_onset_mep_start = auc_start
                    csp_onset_mep_end = auc_end
                    
                    # Init outputs
                    csp_sd_thresh, csp_mcd_thresh = np.nan, np.nan
                    csp_end_sd, csp_end_mcd = None, None
                    # Durations (dict for easier handling)
                    durs = {k: np.nan for k in ['sd_pulse', 'sd_mep_on', 'sd_mep_off', 'mcd_pulse', 'mcd_mep_on', 'mcd_mep_off']}

                    if not short_pre_df.empty:
                        # Thresholds
                        base_rect = short_pre_df[emg_column_name].abs()
                        csp_sd_thresh = csp_sd_multiplier * base_rect.std()
                        
                        mcd_val = base_rect.diff().abs().mean()
                        csp_mcd_thresh = base_rect.mean() - (csp_mcd_multiplier * mcd_val)
                        if csp_mcd_thresh < 0: csp_mcd_thresh = 0
                        
                        # Search
                        search_start = pulse_time + (mep_window_ms[1]/1000.0) + 0.05 # +50ms buffer
                        search_end = pulse_time + (csp_search_window_ms/1000.0)
                        
                        search_df = section_df[(section_df[time_column_name] > search_start) & (section_df[time_column_name] <= search_end)].copy()
                        
                        if not search_df.empty:
                            smoothed = search_df[rectified_emg_col].rolling(window=smoothing_samples_needed).mean()
                            
                            # -- SD Method --
                            above_sd = (smoothed > csp_sd_thresh)
                            cons_sd = above_sd.rolling(window=duration_samples_needed).sum()
                            # Find first index where we have been above threshold for 'duration' samples
                            end_idx_sd = cons_sd[cons_sd >= duration_samples_needed].first_valid_index()
                            
                            if end_idx_sd:
                                # Backtrack to find exact start of return
                                loc = search_df.index.get_loc(end_idx_sd)
                                start_loc = max(0, loc - duration_samples_needed + 1)
                                csp_end_sd = search_df.iloc[start_loc][time_column_name]
                                
                                durs['sd_pulse'] = (csp_end_sd - csp_onset_pulse) * 1000
                                if csp_onset_mep_start: durs['sd_mep_on'] = (csp_end_sd - csp_onset_mep_start) * 1000
                                if csp_onset_mep_end: durs['sd_mep_off'] = (csp_end_sd - csp_onset_mep_end) * 1000
                            
                            # -- MCD Method --
                            above_mcd = (smoothed > csp_mcd_thresh)
                            cons_mcd = above_mcd.rolling(window=duration_samples_needed).sum()
                            end_idx_mcd = cons_mcd[cons_mcd >= duration_samples_needed].first_valid_index()
                            
                            if end_idx_mcd:
                                loc = search_df.index.get_loc(end_idx_mcd)
                                start_loc = max(0, loc - duration_samples_needed + 1)
                                csp_end_mcd = search_df.iloc[start_loc][time_column_name]
                                
                                durs['mcd_pulse'] = (csp_end_mcd - csp_onset_pulse) * 1000
                                if csp_onset_mep_start: durs['mcd_mep_on'] = (csp_end_mcd - csp_onset_mep_start) * 1000
                                if csp_onset_mep_end: durs['mcd_mep_off'] = (csp_end_mcd - csp_onset_mep_end) * 1000

                    # E. Normalization
                    mep_amp_norm, mep_auc_norm = np.nan, np.nan
                    if m_max_amp_used and m_max_amp_used > 0:
                        mep_amp_norm = (mep_amplitude / m_max_amp_used) * 100
                    if m_max_auc_used and m_max_auc_used > 0:
                        mep_auc_norm = (mep_auc / m_max_auc_used) * 100
                        
                    # F. Append Results
                    res = {
                        'filename': file_path.name, 'section': section_num, 'pulse_time_s': pulse_time,
                        'mep_amplitude': mep_amplitude, 'mep_auc': mep_auc,
                        'mep_amplitude_percent_mmax': mep_amp_norm, 'mep_auc_percent_mmax': mep_auc_norm,
                        'mmax_amp_used': m_max_amp_used, 'mmax_auc_used': m_max_auc_used,
                        'short_pre_stim_emg_rms': short_rms, 'long_pre_stim_emg_rms': long_rms,
                        'csp_threshold_sd': csp_sd_thresh, 'csp_threshold_mcd': csp_mcd_thresh,
                        'csp_end_time_sd': csp_end_sd, 'csp_end_time_mcd': csp_end_mcd,
                        'csp_duration_sd_pulse': durs['sd_pulse'], 'csp_duration_sd_mep_onset': durs['sd_mep_on'], 'csp_duration_sd_mep_offset': durs['sd_mep_off'],
                        'csp_duration_mcd_pulse': durs['mcd_pulse'], 'csp_duration_mcd_mep_onset': durs['mcd_mep_on'], 'csp_duration_mcd_mep_offset': durs['mcd_mep_off'],
                    }
                    if local_force_column:
                        res.update({'short_pre_stim_mean_force': short_mean_force, 
                                    'long_pre_stim_mean_force': long_mean_force, 'long_pre_stim_force_cv': long_force_cv})
                    all_results.append(res)
                    
                    # G. Plotting
                    # Raw
                    if plot_rectified_emg:
                         ax_emg.plot(plot_window_df[time_column_name], plot_window_df[rectified_emg_col], color='lightgrey', lw=0.5)
                    ax_emg.plot(plot_window_df[time_column_name], plot_window_df[emg_column_name], color='black', lw=0.5)
                    
                    # Markers
                    ax_emg.axvline(pulse_time, color='red', linestyle='--', label='TMS')
                    ax_emg.axhline(csp_sd_thresh, color='blue', linestyle=':', lw=1, label='SD Thresh')
                    ax_emg.axhline(csp_mcd_thresh, color='cyan', linestyle=':', lw=1, label='MCD Thresh')
                    
                    if csp_end_sd: ax_emg.axvline(csp_end_sd, color='orange', linestyle='--', label='SD End')
                    if csp_end_mcd: ax_emg.axvline(csp_end_mcd, color='purple', linestyle='--', label='MCD End')
                    if mep_max_time: 
                        ax_emg.plot(mep_max_time, mep_max_val, 'bo', ms=4)
                        ax_emg.plot(mep_min_time, mep_min_val, 'ro', ms=4)
                    
                    ax_emg.set_title(f"TMS at {pulse_time:.3f}s | MEP: {mep_amplitude:.2f} | Norm: {mep_amp_norm:.1f}%", fontsize=10)

                # --- Finalize Axes ---
                ax_emg.set_ylabel(emg_column_name)
                ax_emg.legend(loc='upper right', fontsize=6)
                
                if local_force_column:
                    ax_force = axes[ax_idx+1]
                    ax_force.plot(plot_window_df[time_column_name], plot_window_df[local_force_column], color='green', lw=0.5)
                    ax_force.set_ylabel('Force')
                    ax_force.grid(True, alpha=0.3)
                
            # --- Save Figure ---
            fig.tight_layout()
            plot_name = f"{filename_no_ext}_Section_{section_num}.jpg"
            fig.savefig(os.path.join(output_plot_path, plot_name), dpi=150)
            plt.close(fig)
            print(f"  - Saved plot: {plot_name}")

    print("\n✅ Analysis complete.")
    return pd.DataFrame(all_results)
