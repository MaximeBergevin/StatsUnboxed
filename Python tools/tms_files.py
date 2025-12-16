import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import numpy as np
import warnings
from scipy import signal

# =========================================================
# 1. PARSER (Standard LabChart Text File Parser)
# =========================================================
def read_labchart_txt(file_path: str, delimiter: str = '\t') -> pd.DataFrame:
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    data_start_index = 0
    for i, line in enumerate(lines):
        stripped_line = line.strip()
        if stripped_line and (stripped_line[0].isdigit() or stripped_line.startswith((',', '-'))):
            data_start_index = i
            break
            
    # Parse Header
    channel_names = []
    metadata_lines = lines[:data_start_index]
    header_line = next((l for l in metadata_lines if 'ChannelTitle=' in l), None)
    
    if header_line:
        header_text = header_line.split('=', 1)[1].strip()
        channel_names = [name.strip() for name in header_text.split(delimiter) if name.strip()]
    
    if not channel_names:
        if data_start_index >= len(lines): return pd.DataFrame()
        first_data_row_parts = lines[data_start_index].strip().split(delimiter)
        num_parts = len(first_data_row_parts)
        if num_parts <= 1: return pd.DataFrame()
        
        if num_parts > 1 and first_data_row_parts[-1].startswith('#*'):
            num_data_cols = num_parts - 2 
        else:
            num_data_cols = num_parts - 1 
        channel_names = [f'Channel_{i+1}' for i in range(num_data_cols)]

    final_col_names = ['Time'] + channel_names + ['Comment']
    num_expected_cols = len(final_col_names)

    data_rows = []
    for line in lines[data_start_index:]:
        line = line.strip()
        if not line: continue
        parts = line.split(delimiter)
        row_len = len(parts)
        
        if row_len == num_expected_cols:
            data_rows.append(parts)
        elif row_len == num_expected_cols - 1:
            parts.append(None)
            data_rows.append(parts)
        elif row_len > num_expected_cols:
            data_part = parts[:num_expected_cols - 1]
            comment_part = " ".join(parts[num_expected_cols - 1:])
            data_rows.append(data_part + [comment_part])
        elif row_len < num_expected_cols - 1:
            parts.extend([None] * (num_expected_cols - row_len))
            data_rows.append(parts)

    if not data_rows: return pd.DataFrame()
    df = pd.DataFrame(data_rows, columns=final_col_names)
    for col in df.columns:
        if col != 'Comment':
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(',', '.', regex=False), errors='coerce')
    df['Comment'] = df['Comment'].astype('string')
    df.dropna(subset=['Time'], inplace=True)
    return df.reset_index(drop=True)

# =========================================================
# 2. ANALYSIS WORKER (MCD ONLY + Full Metrics)
# =========================================================
def _analyze_single_file(
    file_path: Path,
    output_plot_folder: Path,
    participant_id: str,
    condition: str,
    calculate_csp: bool,
    emg_column_name: str,
    force_column_name: str = None,
    m_wave_trigger_text: str = None,
    tms_trigger_text: str = '#* TMS Pulse',
    time_column_name: str = 'Time',
    # Tuning Params
    pre_stim_baseline_window_ms=(-100, -10),
    mep_window_ms=(5.0, 50.0),
    csp_mcd_multiplier=2.66,       # Default MCD multiplier
    csp_search_window_ms=400.0,
    csp_smoothing_window_ms=10.0,
    csp_min_return_duration_ms=5.0,
    mep_auc_window_ms=(15.0, 100.0),
    mep_onset_sd_thresh=3.0,
    plot_rectified_emg=False
) -> list:
    
    os.makedirs(output_plot_folder, exist_ok=True)
    file_results = []

    try:
        df = read_labchart_txt(str(file_path))
        if df.empty or time_column_name not in df.columns or emg_column_name not in df.columns:
            return []
        local_force_col = force_column_name if (force_column_name and force_column_name in df.columns) else None
    except Exception as e:
        print(f"    ❌ Error reading {file_path.name}: {e}")
        return []

    # --- Section Splitting ---
    time_diffs = df[time_column_name].diff()
    median_step = time_diffs.median()
    if pd.isna(median_step) or median_step <= 0:
        section_starts = [df.index[0]]
    else:
        threshold = median_step * 100
        discontinuities = list(time_diffs[(time_diffs < 0) | (time_diffs > threshold)].index)
        section_starts = sorted(list(set([df.index[0]] + discontinuities)))
    
    sections = []
    for start, end in zip(section_starts, section_starts[1:] + [None]):
        sections.append(df.loc[start:end-1] if end else df.loc[start:])

    # --- Pre-calculate M-waves (Amplitude & AUC) ---
    section_m_waves = {}
    valid_m_sections = []
    if m_wave_trigger_text:
        for s_idx, s_df in enumerate(sections):
            m_indices = s_df[s_df['Comment'].astype(str).str.contains(m_wave_trigger_text, na=False)].index
            if not m_indices.empty:
                pulse_time = s_df.loc[m_indices[-1], time_column_name]
                mw_df = s_df[(s_df[time_column_name] >= pulse_time) & (s_df[time_column_name] <= pulse_time + 0.050)]
                
                if not mw_df.empty:
                    # Amplitude
                    amp = mw_df[emg_column_name].max() - mw_df[emg_column_name].min()
                    
                    # AUC
                    auc = 0
                    base_df = s_df[(s_df[time_column_name] >= pulse_time - 0.100) & (s_df[time_column_name] < pulse_time)]
                    if not base_df.empty:
                        base_mean = base_df[emg_column_name].abs().mean()
                        rect_sig = mw_df[emg_column_name].abs()
                        corr_sig = rect_sig - base_mean
                        corr_sig[corr_sig < 0] = 0 
                        auc = np.trapz(y=corr_sig, x=mw_df[time_column_name])
                    
                    section_m_waves[s_idx] = {'amplitude': amp, 'auc': auc}
                    valid_m_sections.append(s_idx)

    # --- Main Analysis Loop ---
    for section_num, section_df in enumerate(sections, 1):
        s_idx = section_num - 1
        if section_df.empty: continue
        
        section_df = section_df.copy()
        rectified_emg_col = f"{emg_column_name}_rectified"
        section_df[rectified_emg_col] = section_df[emg_column_name].abs()

        # Calculate sampling metrics
        dt = section_df[time_column_name].diff().median()
        if pd.isna(dt) or dt <= 0: dt = 0.001 
        
        smoothing_samples = int(round((csp_smoothing_window_ms / 1000.0) / dt)) or 1
        duration_samples = int(round((csp_min_return_duration_ms / 1000.0) / dt)) or 1

        # Collect Stimuli
        all_stimuli = []
        if m_wave_trigger_text:
            m_indices = section_df[section_df['Comment'].astype(str).str.contains(m_wave_trigger_text, na=False)].index
            for index in m_indices:
                all_stimuli.append({'time': section_df.loc[index, time_column_name], 'type': 'M-Wave'})
        
        tms_indices = section_df[section_df['Comment'].astype(str).str.contains(tms_trigger_text, na=False)].index
        for index in tms_indices:
            all_stimuli.append({'time': section_df.loc[index, time_column_name], 'type': 'TMS'})
        
        if not all_stimuli: continue
        all_stimuli.sort(key=lambda x: x['time'])

        # Get M-Max for this section
        m_max_amp_used, m_max_auc_used = None, None
        if s_idx in section_m_waves:
            m_max_amp_used = section_m_waves[s_idx]['amplitude']
            m_max_auc_used = section_m_waves[s_idx]['auc']
        elif valid_m_sections:
            closest_s_idx = min(valid_m_sections, key=lambda x: abs(x - s_idx))
            m_max_amp_used = section_m_waves[closest_s_idx]['amplitude']
            m_max_auc_used = section_m_waves[closest_s_idx]['auc']

        # Plot Initialization
        plots_per_pulse = 2 if local_force_col else 1
        fig_height = max(4, (3 * len(all_stimuli)) * plots_per_pulse)
        fig, axes = plt.subplots(nrows=len(all_stimuli) * plots_per_pulse, ncols=1, figsize=(12, fig_height), squeeze=False)
        axes = axes.flatten()

        for i, stim in enumerate(all_stimuli):
            pulse_time = stim['time']
            stim_type = stim['type']
            ax_emg = axes[i * plots_per_pulse]
            
            # Plot window
            plot_window = section_df[(section_df[time_column_name] >= pulse_time - 0.2) & (section_df[time_column_name] <= pulse_time + 0.5)]

            # ---------------------------
            # 1. TYPE: M-WAVE
            # ---------------------------
            if stim_type == 'M-Wave':
                ax_emg.plot(plot_window[time_column_name], plot_window[emg_column_name], color='black', lw=0.5)
                ax_emg.axvline(pulse_time, color='blue', linestyle='--', label='M-Wave Stim')
                title = f"M-Wave (Section {section_num})"
                
                # Check if we have calculated stats for this specific pulse to display
                # (We already calculated m-max for the section, let's just display that)
                if m_max_amp_used:
                    title += f" | Amp: {m_max_amp_used:.2f}mV | AUC: {m_max_auc_used:.2f}"
                ax_emg.set_title(title, fontsize=10)

            # ---------------------------
            # 2. TYPE: TMS
            # ---------------------------
            elif stim_type == 'TMS':
                # --- A. MEP AMPLITUDE ---
                mep_start = pulse_time + (mep_window_ms[0] / 1000.0)
                mep_end = pulse_time + (mep_window_ms[1] / 1000.0)
                mep_df = section_df[(section_df[time_column_name] >= mep_start) & (section_df[time_column_name] <= mep_end)]
                
                mep_amplitude = 0
                if not mep_df.empty:
                    mep_amplitude = mep_df[emg_column_name].max() - mep_df[emg_column_name].min()

                # --- B. BASELINES ---
                short_pre_start = pulse_time + (pre_stim_baseline_window_ms[0]/1000.0)
                short_pre_end = pulse_time + (pre_stim_baseline_window_ms[1]/1000.0)
                short_pre_df = section_df[(section_df[time_column_name] >= short_pre_start) & (section_df[time_column_name] < short_pre_end)]
                
                short_rms = np.nan
                if not short_pre_df.empty:
                    short_rms = np.sqrt(np.mean(short_pre_df[emg_column_name]**2))

                # --- C. MEP AUC (Complex Logic) ---
                mep_auc = 0
                auc_start, auc_end = None, None
                if not short_pre_df.empty:
                    base_mean = short_pre_df[emg_column_name].mean()
                    base_sd = short_pre_df[emg_column_name].std()
                    upper = base_mean + mep_onset_sd_thresh * base_sd
                    lower = base_mean - mep_onset_sd_thresh * base_sd
                    
                    # Onset (20-50ms)
                    onset_search = section_df[(section_df[time_column_name] >= pulse_time + 0.020) & (section_df[time_column_name] <= pulse_time + 0.050)]
                    onset_pts = onset_search[(onset_search[emg_column_name] > upper) | (onset_search[emg_column_name] < lower)]
                    auc_start = onset_pts.iloc[0][time_column_name] if not onset_pts.empty else None
                    
                    # Offset (50-150ms)
                    if auc_start:
                         offset_search = section_df[(section_df[time_column_name] >= pulse_time + 0.050) & (section_df[time_column_name] <= pulse_time + 0.150)]
                         offset_pts = offset_search[(offset_search[emg_column_name] < upper) & (offset_search[emg_column_name] > lower)]
                         auc_end = offset_pts.iloc[0][time_column_name] if not offset_pts.empty else None

                    # Fallback
                    if auc_start is None or auc_end is None:
                        auc_start = pulse_time + (mep_auc_window_ms[0]/1000.0)
                        auc_end = pulse_time + (mep_auc_window_ms[1]/1000.0)
                        
                    auc_df = section_df[(section_df[time_column_name] >= auc_start) & (section_df[time_column_name] <= auc_end)]
                    if not auc_df.empty:
                         rect_base = short_pre_df[emg_column_name].abs().mean()
                         rect_mep = auc_df[emg_column_name].abs()
                         corr_mep = rect_mep - rect_base
                         corr_mep[corr_mep < 0] = 0
                         mep_auc = np.trapz(y=corr_mep, x=auc_df[time_column_name])

                # --- D. CSP Calculations (MCD ONLY) ---
                csp_end_mcd = None
                csp_duration_mcd = np.nan
                csp_mcd_thresh = np.nan

                if calculate_csp and not short_pre_df.empty:
                    base_rect = short_pre_df[emg_column_name].abs()
                    
                    # MCD Calculation
                    mcd_val = base_rect.diff().abs().mean()
                    csp_mcd_thresh = base_rect.mean() - (csp_mcd_multiplier * mcd_val)
                    if csp_mcd_thresh < 0: csp_mcd_thresh = 0
                    
                    # Search Window
                    search_start = pulse_time + (mep_window_ms[1]/1000.0) + 0.05
                    search_end = pulse_time + (csp_search_window_ms/1000.0)
                    
                    search_df = section_df[(section_df[time_column_name] > search_start) & (section_df[time_column_name] <= search_end)].copy()
                    
                    if not search_df.empty:
                        smoothed = search_df[rectified_emg_col].rolling(window=smoothing_samples).mean()
                        
                        # -- MCD Logic --
                        above_mcd = (smoothed > csp_mcd_thresh)
                        # Check for CONSECUTIVE samples above threshold
                        cons_mcd = above_mcd.rolling(window=duration_samples).sum()
                        
                        # The first index where the sum == window size is the END of the return to baseline
                        end_idx_mcd = cons_mcd[cons_mcd >= duration_samples].first_valid_index()
                        
                        if end_idx_mcd:
                            # Backtrack to the start of that window
                            loc = search_df.index.get_loc(end_idx_mcd)
                            start_loc = max(0, loc - duration_samples + 1)
                            csp_end_mcd = search_df.iloc[start_loc][time_column_name]
                            csp_duration_mcd = (csp_end_mcd - pulse_time) * 1000

                # --- E. Normalization ---
                mep_amp_norm, mep_auc_norm = np.nan, np.nan
                if m_max_amp_used and m_max_amp_used > 0:
                    mep_amp_norm = (mep_amplitude / m_max_amp_used) * 100
                if m_max_auc_used and m_max_auc_used > 0:
                    mep_auc_norm = (mep_auc / m_max_auc_used) * 100

                # --- F. Append Results ---
                res = {
                    'Participant': participant_id, 
                    'Condition': condition, 
                    'Filename': file_path.name,
                    'Section': section_num, 
                    'Pulse_Time': pulse_time,
                    'MEP_Amplitude_mV': mep_amplitude, 
                    'MEP_AUC': mep_auc,
                    'MEP_Amplitude_Norm_Percent': mep_amp_norm, 
                    'MEP_AUC_Norm_Percent': mep_auc_norm,
                    'M_Max_Amp_Used': m_max_amp_used,
                    'M_Max_AUC_Used': m_max_auc_used,
                    'Baseline_RMS': short_rms,
                    'CSP_End_Time_MCD': csp_end_mcd,
                    'CSP_Duration_MCD_ms': csp_duration_mcd,
                    'CSP_Threshold_MCD': csp_mcd_thresh
                }
                if local_force_col:
                    res['Mean_Force'] = plot_window[local_force_col].mean()
                
                file_results.append(res)

                # --- G. Plotting ---
                if plot_rectified_emg:
                     ax_emg.plot(plot_window[time_column_name], plot_window[rectified_emg_col], color='lightgrey', lw=0.5)
                ax_emg.plot(plot_window[time_column_name], plot_window[emg_column_name], color='black', lw=0.5)
                
                ax_emg.axvline(pulse_time, color='red', linestyle='--', label='TMS')
                
                if calculate_csp and not np.isnan(csp_mcd_thresh):
                    ax_emg.axhline(csp_mcd_thresh, color='cyan', linestyle=':', lw=1.5, label='MCD Thresh')
                    if csp_end_mcd: 
                        ax_emg.axvline(csp_end_mcd, color='purple', linestyle='--', lw=1.5, label='MCD End')

                title = f"TMS ({condition}) | MEP: {mep_amplitude:.2f}mV"
                if calculate_csp: 
                    title += f" | CSP (MCD): {csp_duration_mcd:.1f}ms"
                ax_emg.set_title(title, fontsize=10)
                ax_emg.legend(loc='upper right', fontsize=6)

            # Force Plot
            if local_force_col:
                ax_force = axes[i*plots_per_pulse+1]
                ax_force.plot(plot_window[time_column_name], plot_window[local_force_col], color='green', lw=0.5)
                ax_force.set_ylabel('Force')
        
        fig.tight_layout()
        plot_name = f"{file_path.stem}_Sec{section_num}.jpg"
        fig.savefig(output_plot_folder / plot_name, dpi=100)
        plt.close(fig)

    return file_results

# =========================================================
# 3. BATCH MANAGER (Auto-iterates Active/Rest)
# =========================================================
def batch_process_study(
    root_folder: str,
    output_folder: str,
    **kwargs
) -> pd.DataFrame:
    
    root_path = Path(root_folder)
    out_path = Path(output_folder)
    
    all_study_data = []
    conditions = ['Active_session', 'Rest_sessions']
    
    print(f"📂 Starting Batch Analysis in: {root_path}")

    for condition in conditions:
        cond_path = root_path / condition
        if not cond_path.exists(): continue
        
        do_csp = (condition == 'Active_session')
        
        participant_folders = [p for p in cond_path.iterdir() if p.is_dir()]
        
        for p_folder in sorted(participant_folders):
            p_id = p_folder.name
            print(f"  👉 Processing {condition} / {p_id} ...")
            
            p_plot_folder = out_path / condition / p_id
            participant_data = []

            txt_files = list(p_folder.glob('*.txt'))
            for txt_file in txt_files:
                results = _analyze_single_file(
                    file_path=txt_file,
                    output_plot_folder=p_plot_folder,
                    participant_id=p_id,
                    condition=condition,
                    calculate_csp=do_csp,
                    **kwargs
                )
                participant_data.extend(results)
                all_study_data.extend(results)
            
            # --- SAVE INDIVIDUAL CSV FOR THIS PARTICIPANT ---
            if participant_data:
                df_p = pd.DataFrame(participant_data)
                csv_name = p_plot_folder / f"{p_id}_{condition}_results.csv"
                df_p.to_csv(csv_name, index=False)
                print(f"     -> Saved CSV: {csv_name.name}")

    print("\n✅ Batch Analysis Complete.")
    return pd.DataFrame(all_study_data)
