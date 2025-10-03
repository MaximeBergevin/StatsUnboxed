import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import numpy as np
import warnings
from scipy import signal

def read_labchart_txt(file_path: str, delimiter: str = '\t') -> pd.DataFrame:
    """
    Reads a LabChart .txt file by manually parsing it line-by-line. This
    function is tailored to handle an unnamed Time column, named data channels, 
    and a final, sparse Comment column. It also automatically removes rows 
    with corrupted or missing time data, including metadata blocks inserted mid-file.

    Args:
        file_path (str): The full path to the .txt file.
        delimiter (str, optional): The column delimiter. Defaults to '\t'.

    Returns:
        pd.DataFrame: A cleaned DataFrame with proper headers and data types.
    """
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # Find the start of the actual data
    # ----------------------------------
    data_start_index = 0
    for i, line in enumerate(lines):
        stripped_line = line.strip()
        if stripped_line and (stripped_line[0].isdigit() or stripped_line.startswith((',', '-'))):
            data_start_index = i
            break
    if data_start_index == 0:
        warnings.warn(f"Could not detect metadata in {file_path}. Assuming data starts on line 1.")

    # Parse the Header to get channel names
    # ---------------------------------------
    channel_names = []
    metadata_lines = lines[:data_start_index]
    header_line = next((l for l in metadata_lines if 'ChannelTitle=' in l), None)
    
    if header_line:
        header_text = header_line.split('=', 1)[1].strip()
        channel_names = [name.strip() for name in header_text.split(delimiter) if name.strip()]
    else:
        warnings.warn("'ChannelTitle=' not found. Header names will be generic.")

    # Process the Data Rows
    # -----------------------
    data_rows = []
    for line in lines[data_start_index:]:
        line = line.strip()
        if not line:
            continue
        
        if '#*' in line:
            parts = line.split(delimiter)
            data_rows.append(parts)
        else:
            parts = line.split(delimiter)
            parts.append(None)
            data_rows.append(parts)

    # Create the DataFrame
    # ----------------------
    if not data_rows:
        return pd.DataFrame()

    if not channel_names:
        channel_names = [f'Channel_{i+1}' for i in range(8)]
        
    final_col_names = ['Time'] + channel_names + ['Comment']
    df = pd.DataFrame(data_rows, columns=final_col_names)

    # Convert data types and clean the data
    # ---------------------------------------
    for col in df.columns:
        if col != 'Comment':
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(',', '.'), errors='coerce')
    df['Comment'] = df['Comment'].astype('string')
    
    df.dropna(subset=['Time'], inplace=True)

    return df.reset_index(drop=True)


def tms_csp(
    folder_path: str,
    emg_column_name: str,
    output_plot_path: str,
    force_column_name: str = None,
    time_column_name: str = 'Time',
    tms_trigger_text: str = '#* TMS Pulse',
    pre_stim_baseline_window_ms: tuple = (-120.0, -20.0),
    long_pre_stim_window_s: float = 2.5,
    mep_window_ms: tuple = (5.0, 50.0),
    csp_threshold_percentile: float = 0.50,
    csp_search_window_ms: float = 400.0,
    mep_auc_window_ms: tuple = (15.0, 100.0),
    mep_onset_sd_thresh: float = 3.0
) -> pd.DataFrame:
    """
    Analyzes LabChart .txt files, calculating MEP, CSP, pre-stimulus RMS (short and long window), 
    and optionally processing a force channel (including mean and CV).
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
            df = read_labchart_txt(file_path)
            local_force_column = force_column_name
            if df.empty or time_column_name not in df.columns or emg_column_name not in df.columns:
                 print(f"  - Error: Required columns not found. Skipping file.")
                 continue
            if local_force_column and local_force_column not in df.columns:
                print(f"  - Warning: Force column '{local_force_column}' not found. Proceeding without force analysis.")
                local_force_column = None
        except Exception as e:
            print(f"  - Error reading {file_path.name}: {e}. Skipping file.")
            continue
        
        # Section Splitting
        # -------------------
        first_valid_index = df.index[0]
        time_diffs = df[time_column_name].diff()
        median_step = time_diffs.median()
        
        if pd.isna(median_step) or median_step <= 0:
            section_start_indices = [first_valid_index]
        else:
            discontinuity_threshold = median_step * 100
            discontinuities = list(time_diffs[(time_diffs < 0) | (time_diffs > discontinuity_threshold)].index)
            section_start_indices = [first_valid_index] + discontinuities
        
        section_start_indices = sorted(list(set(section_start_indices)))
        sections = [df.loc[start:end-1] if end is not None else df.loc[start:] for start, end in zip(section_start_indices, section_start_indices[1:] + [None])]
        
        print(f"  - Detected {len(sections)} section(s) in this file.")
        
        # Main Analysis Loop
        # --------------------
        for section_num, section_df in enumerate(sections, 1):
            if section_df.empty:
                continue
            
            section_df = section_df.copy()
            rectified_emg_col = f"{emg_column_name}_rectified"
            section_df[rectified_emg_col] = section_df[emg_column_name].abs()

            pulse_indices = section_df[section_df['Comment'].astype(str).str.contains(tms_trigger_text, na=False)].index
            pulse_times = section_df.loc[pulse_indices, time_column_name].tolist()

            if not pulse_times:
                print(f"  - No TMS pulses found in Section #{section_num}.")
                continue
                
            print(f"  - Section #{section_num}: Found {len(pulse_times)} TMS pulse(s). Analyzing...")
            
            plots_per_pulse = 2 if local_force_column else 1
            num_rows = len(pulse_times) * plots_per_pulse
            fig_height = (3 * len(pulse_times)) * plots_per_pulse
            fig, axes = plt.subplots(nrows=num_rows, ncols=1, figsize=(12, fig_height), squeeze=False)
            axes = axes.flatten()

            for i, pulse_time in enumerate(pulse_times):
                
                # MEP Amplitude Calculation
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

                # Short Pre-stimulus Baseline (for thresholds)
                # ----------------------------------------------
                short_pre_stim_start = pulse_time + (pre_stim_baseline_window_ms[0] / 1000.0)
                short_pre_stim_end = pulse_time + (pre_stim_baseline_window_ms[1] / 1000.0)
                short_pre_stim_df = section_df[(section_df[time_column_name] >= short_pre_stim_start) & (section_df[time_column_name] < short_pre_stim_end)]
                
                short_pre_stim_emg_rms = float('nan')
                short_pre_stim_mean_force = float('nan')

                if short_pre_stim_df.empty:
                    print(f"    - Warning: No short pre-stimulus data for pulse at {pulse_time:.3f}s.")
                else:
                    short_pre_stim_emg_rms = np.sqrt(np.mean(short_pre_stim_df[emg_column_name]**2))
                    if local_force_column:
                        short_pre_stim_mean_force = short_pre_stim_df[local_force_column].mean()

                # Long Pre-stimulus Baseline (for stability measures)
                # ---------------------------------------------------
                long_pre_stim_start = pulse_time - long_pre_stim_window_s
                long_pre_stim_end = pulse_time
                long_pre_stim_df = section_df[(section_df[time_column_name] >= long_pre_stim_start) & (section_df[time_column_name] < long_pre_stim_end)]
                
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
                
                # Dynamic MEP Window and AUC Calculation
                # ---------------------------------------
                rectified_baseline_emg = short_pre_stim_df[emg_column_name].abs()
                pre_stim_mean_rectified = rectified_baseline_emg.mean()
                baseline_mean = short_pre_stim_df[emg_column_name].mean()
                baseline_sd = short_pre_stim_df[emg_column_name].std()
                upper_thresh = baseline_mean + mep_onset_sd_thresh * baseline_sd
                lower_thresh = baseline_mean - mep_onset_sd_thresh * baseline_sd
                onset_search_start = pulse_time + 0.010
                onset_search_end = pulse_time + (mep_window_ms[1] / 1000.0)
                onset_search_df = section_df[(section_df[time_column_name] >= onset_search_start) & (section_df[time_column_name] <= onset_search_end)]
                onset_point = onset_search_df[(onset_search_df[emg_column_name] > upper_thresh) | (onset_search_df[emg_column_name] < lower_thresh)]
                auc_start_time = onset_point.iloc[0][time_column_name] if not onset_point.empty else None
                offset_search_start = pulse_time + (mep_window_ms[1] / 1000.0)
                offset_search_end = pulse_time + 0.150
                offset_search_df = section_df[(section_df[time_column_name] >= offset_search_start) & (section_df[time_column_name] <= offset_search_end)]
                offset_point = offset_search_df[(offset_search_df[emg_column_name] < upper_thresh) & (offset_search_df[emg_column_name] > lower_thresh)]
                auc_end_time = offset_point.iloc[0][time_column_name] if not offset_point.empty else None
                if auc_start_time is None or auc_end_time is None:
                    auc_start_time = pulse_time + (mep_auc_window_ms[0] / 1000.0)
                    auc_end_time = pulse_time + (mep_auc_window_ms[1] / 1000.0)
                mep_auc_df = section_df[(section_df[time_column_name] >= auc_start_time) & (section_df[time_column_name] <= auc_end_time)].copy()
                mep_auc = 0
                if not mep_auc_df.empty:
                    rectified_mep = mep_auc_df[emg_column_name].abs()
                    corrected_rectified_mep = rectified_mep - pre_stim_mean_rectified
                    corrected_rectified_mep[corrected_rectified_mep < 0] = 0
                    mep_auc = np.trapz(y=corrected_rectified_mep, x=mep_auc_df[time_column_name])
                
                # CSP Detection
                # ----------------
                csp_threshold_value = rectified_baseline_emg.quantile(csp_threshold_percentile)
                csp_search_start = pulse_time + (mep_window_ms[1] / 1000.0) + 0.05 # 50ms after MEP window
                csp_search_end = pulse_time + (csp_search_window_ms / 1000.0)
                csp_search_df = section_df[(section_df[time_column_name] > csp_search_start) & (section_df[time_column_name] <= csp_search_end)]
                csp_end_time, csp_duration_ms = None, float('nan')
                if not csp_search_df.empty:
                    activity_return_df = csp_search_df[csp_search_df[rectified_emg_col] > csp_threshold_value]
                    if not activity_return_df.empty:
                        csp_end_time = activity_return_df.iloc[0][time_column_name]
                        csp_duration_ms = (csp_end_time - pulse_time) * 1000.0
                
                # Store Results
                # ---------------
                results_to_append = {
                    'filename': file_path.name, 'section': section_num, 'pulse_time_s': pulse_time,
                    'mep_amplitude': mep_amplitude, 'mep_auc': mep_auc, 'csp_duration_ms': csp_duration_ms,
                    'short_pre_stim_emg_rms': short_pre_stim_emg_rms,
                    'long_pre_stim_emg_rms': long_pre_stim_emg_rms,
                }
                if local_force_column:
                    results_to_append['short_pre_stim_mean_force'] = short_pre_stim_mean_force
                    results_to_append['long_pre_stim_mean_force'] = long_pre_stim_mean_force
                    results_to_append['long_pre_stim_force_cv'] = long_pre_stim_force_cv
                
                all_results.append(results_to_append)

                # Facetted Plotting Logic
                # -------------------------
                ax_idx = i * plots_per_pulse
                ax_emg = axes[ax_idx]
                plot_window_df = section_df[(section_df[time_column_name] >= pulse_time - 0.2) & (section_df[time_column_name] <= pulse_time + 0.5)]
                ax_emg.plot(plot_window_df[time_column_name], plot_window_df[emg_column_name], color='black', lw=0.5)
                ax_emg.set_ylabel(f"EMG Amp ({emg_column_name})")
                ax_emg.axvline(pulse_time, color='red', linestyle='--', label=f'TMS Pulse')
                ax_emg.axhline(csp_threshold_value, color='blue', linestyle=':', lw=1, label=f'CSP Threshold')
                if not mep_auc_df.empty:
                    ax_emg.fill_between(mep_auc_df[time_column_name], 0, mep_auc_df[emg_column_name], color='green', alpha=0.3, label=f'MEP AUC: {mep_auc:.4f}')
                if mep_max_time is not None:
                    ax_emg.plot(mep_max_time, mep_max_val, 'o', color='blue', markersize=5)
                    ax_emg.plot(mep_min_time, mep_min_val, 'o', color='red', markersize=5, label=f'MEP Amp: {mep_amplitude:.2f}')
                if csp_end_time:
                    ax_emg.axvspan(pulse_time, csp_end_time, color='gray', alpha=0.3, label=f'CSP: {csp_duration_ms:.1f} ms')
                    ax_emg.axvline(csp_end_time, color='orange', linestyle='--', label=f'CSP End')
                ax_emg.set_title(f"Pulse at {pulse_time:.3f} s", fontsize=10)
                ax_emg.legend(loc='upper right', fontsize=8)
                ax_emg.grid(True, linestyle='--', alpha=0.6)
                if local_force_column:
                    ax_force = axes[ax_idx + 1]
                    ax_emg.get_shared_x_axes().joined(ax_emg, ax_force)
                    ax_force.plot(plot_window_df[time_column_name], plot_window_df[local_force_column], color='green', lw=0.5, linestyle='--')
                    ax_force.set_ylabel(f"Force ({local_force_column})")
                    ax_force.set_xlabel("Time (s)")
                    ax_force.grid(True, linestyle='--', alpha=0.6)
                    ax_emg.tick_params(axis='x', labelbottom=False)
                else:
                    ax_emg.set_xlabel("Time (s)")
            
            # Save Figure
            # -------------
            fig.suptitle(f'Analysis for: {file_path.name} - Section #{section_num}', fontsize=16, y=1.0)
            fig.tight_layout(rect=[0, 0, 1, 0.98], h_pad=3.0)
            plot_filename = f"{filename_no_ext}_Section_{section_num}.jpg"
            fig.savefig(os.path.join(output_plot_path, plot_filename), dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"  - Plot saved to {plot_filename}")

    print("\n✅ Analysis complete.")
    return pd.DataFrame(all_results)











