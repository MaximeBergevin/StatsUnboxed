import os
import tms_tools 
import pandas as pd

# This ensures the script runs only when executed directly
if __name__ == "__main__":
    
    # =================================================================
    # 1. SETUP: Define paths and analysis parameters
    # =================================================================
    
    base_folder = os.getcwd()
    data_folder = os.path.join(base_folder, 'active_MEP')
    output_folder = os.path.join(base_folder, 'plots')
    
    os.makedirs(output_folder, exist_ok=True)
    
    # --- Static Parameters ---
    emg_channel_name = 'Flex EMG'
    force_channel_name = '%MVC'
    m_wave_id_text = '#* Digitimer'

    print("--- Starting TMS Data Analysis ---")
    print(f"Data folder: {data_folder}")
    print(f"Output folder: {output_folder}")
    print(f"EMG Channel: '{emg_channel_name}'")
    
    # =================================================================
    # 2. DEFINE ANALYSIS CONFIGURATIONS
    # =================================================================

    # --- Config 1: MCD Method (Garvey et al. 2001) ---
    # This is robust to outliers and based on the paper.
    config_mcd = {
        "csp_method": "mcd",
        "csp_mcd_multiplier": 2.66, # 2.66 = 99.76% (3-SD) limit
        "csp_smoothing_window_ms": 10.0,
        "csp_min_return_duration_ms": 100.0,
        "pre_stim_baseline_window_ms": (-100, 0),
        "output_suffix": "mcd_2-66"
    }

    # --- Config 2: Median Method (Robust) ---
    # Use a low percentile (10-30%) to set a sensitive "floor"
    config_median = {
        "csp_method": "median",
        "csp_threshold_percentile": 0.10, # 10th percentile
        "csp_smoothing_window_ms": 10.0,
        "csp_min_return_duration_ms": 1.0,
        "pre_stim_baseline_window_ms": (-100, 0),
        "output_suffix": "median_10pct"
    }
    
    # --- Config 3: SD Method (Outlier-Sensitive) ---
    # This will fail if there are artifacts.
    # Uses the (-200, -100) window to avoid pre-pulse artifacts.
    config_sd = {
        "csp_method": "sd",
        "csp_sd_multiplier": 2.0,
        "csp_smoothing_window_ms": 10.0,
        "csp_min_return_duration_ms": 25.0, # Requires a longer duration
        "pre_stim_baseline_window_ms": (-200, -100), # Artifact fix
        "output_suffix": "sd_2-0_artifact_fix"
    }

    # =================================================================
    # 3. EXECUTION: Select a config and run
    # =================================================================
    
    # --- CHOOSE YOUR ANALYSIS ---
    # Just change this one line to test a different method
    active_config = config_sd
    # active_config = config_median
    # active_config = config_sd
    # ----------------------------

    print(f"\nRunning analysis with '{active_config['csp_method']}' method...")

    # Base parameters (these rarely change)
    base_params = {
        "folder_path": data_folder,
        "output_plot_path": output_folder,
        "emg_column_name": emg_channel_name,
        "force_column_name": force_channel_name,
        "m_wave_trigger_text": m_wave_id_text,
        "mep_window_ms": (5.0, 50.0),
        "mep_auc_window_ms": (15.0, 100.0),
        "csp_search_window_ms": 400.0,
        "plot_rectified_emg": True
    }

    # Combine the base parameters with the active config
    # This uses all key-value pairs from both dictionaries
    analysis_params = {**base_params, **active_config}
    
    # Get the output file suffix and remove it from the dict
    # (tms_csp doesn't accept 'output_suffix' as an argument)
    output_suffix = analysis_params.pop("output_suffix")

    # Run the analysis by unpacking the parameters dictionary
    results_df = tms_tools.tms_csp(**analysis_params)
    
    # =================================================================
    # 4. OUTPUT: Display and save the results
    # =================================================================
    
    if results_df is not None and not results_df.empty:
        print(f"\n--- Analysis Complete: First 5 results ({results_df.shape[0]} total) ---")
        print(results_df.head())

        # Use the suffix to create unique filenames
        csv_path = os.path.join(output_folder, f'tms_analysis_results_{output_suffix}.csv')
        excel_path = os.path.join(output_folder, f'tms_analysis_results_{output_suffix}.xlsx')

        # Save to CSV
        try:
            results_df.to_csv(csv_path, index=False)
            print(f"\n✅ Results table saved to: {csv_path}")
        except Exception as e:
            print(f"\n❌ Could not save CSV file: {e}")
        
        # Save to Excel
        try:
            results_df.to_excel(excel_path, index=False)
            print(f"✅ Results table also saved to: {excel_path}")
        except Exception as e:
            print(f"\n❌ Could not save Excel file: {e}. (Is openpyxl installed?)")
            
    else:
        print("\n--- Analysis finished, but no results were generated. ---")

        print("Please check the console for warnings about missing files or columns.")
