import pandas as pd
import os
import glob

def merge_quant_files(study_number, run_numbers, output_dir="Output"):
    quant_dfs = []  # Corrected the variable name
    
    for run_number in run_numbers:
        # Find directories matching the pattern
        run_dir_pattern = os.path.join(output_dir, study_number, f"{run_number}_pass", "quant", f"{run_number}_pass")
        matching_dirs = glob.glob(run_dir_pattern)
        
        if not matching_dirs:
            print(f"No matching directory found for pattern: {run_dir_pattern}")
            raise FileNotFoundError(f"No matching directory found for pattern: {run_dir_pattern}")
        
        for run_dir in matching_dirs:
            quant_file_path = os.path.join(run_dir, "quant.sf")
            if os.path.exists(quant_file_path):
                df = pd.read_csv(quant_file_path, sep='\t', usecols=['Name', 'TPM'])
                df.rename(columns={'TPM': run_number}, inplace=True)
                quant_dfs.append(df)
                print(f"File {quant_file_path} found and added. DataFrame shape: {df.shape}")
            else:
                print(f"File {quant_file_path} not found.")
                raise FileNotFoundError(f"File {quant_file_path} not found.")
    
    if quant_dfs:
        merged_df = quant_dfs[0]
        for df in quant_dfs[1:]:
            merged_df = pd.merge(merged_df, df, on='Name')
        
        merged_file_path = os.path.join(output_dir, study_number, "merged_quant.sf")
        os.makedirs(os.path.dirname(merged_file_path), exist_ok=True)
        merged_df.to_csv(merged_file_path, sep='\t', index=False)
        print(f"Merged file created at {merged_file_path}. Merged DataFrame shape: {merged_df.shape}")
        return merged_file_path
    else:
        raise ValueError("No quant.sf files found to merge.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Merge quant.sf files for a study.")
    parser.add_argument("study_number", type=str, help="The study number.")
    parser.add_argument("run_numbers", nargs='+', type=str, help="The run numbers for the study.")
    parser.add_argument("--output_dir", type=str, default="Output", help="The output directory where the merged file will be saved.")
    
    args = parser.parse_args()
    
    try:
        merge_quant_files(args.study_number, args.run_numbers, args.output_dir)
    except Exception as e:
        print(f"Error: {e}")
