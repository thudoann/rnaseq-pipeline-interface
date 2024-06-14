# Import required libraries
import pandas as pd
import subprocess
import os
import shutil
import gzip
from collections import defaultdict

# Read the CSV file into a DataFrame
df = pd.read_csv('SraRunTable.csv', delimiter=',')

# Create the Study folder
study_folder = "Study"
os.makedirs(study_folder, exist_ok=True)

# Extract unique SRA_study names
unique_sra_studies = df['SRA_study'].unique()

for sra_study in unique_sra_studies:
    # Convert sra_study to string
    sra_study = str(sra_study)
    
    # Create folder for each SRA_study name inside Study folder
    sra_study_folder = os.path.join(study_folder, sra_study)
    os.makedirs(sra_study_folder, exist_ok=True)
    
    # Directory to store SRA files for this SRA_study
    sra_directory = os.path.join(sra_study_folder, "SRA")
    os.makedirs(sra_directory, exist_ok=True)    
    

    # Directory to store fastq files
    fastq_directory = os.path.join(sra_study_folder, "fastq")
    if not os.path.exists(fastq_directory):  # Check if the directory already exists
        os.makedirs(fastq_directory)
    
    # Directory to store trimmed fastq files
    trimmed_fastq_directory = os.path.join(fastq_directory, "trimmed_fastq")
    os.makedirs(trimmed_fastq_directory, exist_ok=True)

    # Download SRA files
    sra_numbers = df[df['SRA_study'] == sra_study]['Run'].tolist()
    for sra_id in sra_numbers:
        print("Currently downloading:", sra_id)
        prefetch_command = ["prefetch", sra_id]
        print("The command used was:", prefetch_command)
        subprocess.run(prefetch_command, cwd=sra_directory)

    
    # Extract fastq files
    for sra_id in sra_numbers:
        sra_id_folder = os.path.join(sra_study_folder, "SRA", sra_id)
        sra_file_path = os.path.join(sra_id_folder, f"{sra_id}.sra")
        if os.path.exists(sra_file_path):
            print("Generating fastq for:", sra_id)
            fastq_dump_command = [
                "fastq-dump",
                "--outdir", fastq_directory,
                "--gzip",
                "--skip-technical",
                "--readids",
                "--read-filter", "pass",
                "--dumpbase",
                "--split-3",
                "--clip",
                sra_file_path
            ]
            print("The command used was:", fastq_dump_command)
            subprocess.run(fastq_dump_command)
        else:
            print("SRA file not found:", sra_file_path)

    # Unzip .fastq.gz files
    for fastq_file in os.listdir(fastq_directory):
        if fastq_file.endswith(".fastq.gz"):
            fastq_file_path = os.path.join(fastq_directory, fastq_file)
            unzipped_file_path = os.path.join(fastq_directory, fastq_file[:-3])  # Remove .gz extension
            with gzip.open(fastq_file_path, 'rb') as f_in:
                with open(unzipped_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(fastq_file_path)  # Remove the original .fastq.gz file


    # Identify the FASTQ file(s)
    fastq_files = [file for file in os.listdir(fastq_directory) if file.endswith(".fastq")]

    # Group FASTQ files by run number
    fastq_files_by_run = defaultdict(list)
    for fastq_file in fastq_files:
        run_number = fastq_file.split('_')[0]  # Assuming the run number is at the beginning of the file name
        fastq_files_by_run[run_number].append(fastq_file)


    # Process each group of FASTQ files
    for run_number, files_in_run in fastq_files_by_run.items():
        # Create a new folder for each group of run numbers
        run_number_folder = os.path.join(trimmed_fastq_directory, run_number)
        os.makedirs(run_number_folder, exist_ok=True)
        # Directory to store SRA files for this SRA_study
        sra_final_dir = os.path.join(sra_study_folder, run_number)
        os.makedirs(sra_final_dir, exist_ok=True)   
 
        # Process based on the number of FASTQ files in the group
        if len(files_in_run) == 1:
            # Run FastQC on each FASTQ file
            for fastq_file in files_in_run:
                fastqc_command = ["fastqc", os.path.join(fastq_directory, fastq_file)]
                print("Running FastQC on:", fastq_file)
                subprocess.run(fastqc_command)
                html_files = [file for file in os.listdir(fastq_directory) if file.endswith(".html")]
                # Move the HTML file to the sra_study_folder
                html_file_path1 = os.path.join(fastq_directory, html_files[0])
                shutil.move(html_file_path1, os.path.join(sra_final_dir, html_files[0]))


            # Run Trim Galore on each FASTQ file
            for fastq_file in files_in_run:
                trim_galore_command = ["trim_galore", os.path.join(fastq_directory, fastq_file), "--output_dir", run_number_folder]
                print("Running Trim Galore on:", fastq_file)
                subprocess.run(trim_galore_command)
                txt_files = [file for file in os.listdir(run_number_folder) if file.endswith(".txt")]
                # Move the txt file to the sra_study_folder
                txt_file_path1 = os.path.join(run_number_folder, txt_files[0])

                shutil.move(txt_file_path1, os.path.join(sra_final_dir, txt_files[0]))


             # Identify the trimmed FASTQ file(s)
            trimmed_fastq_files = [file for file in os.listdir(run_number_folder) if file.endswith(".fq")]

           

            # Run Salmon indexing
            salmon_index_command = ["salmon", "index", "-t", "./t_indxs/pao1_cnda.fa.gz", "-i", "./t_indxs/pao1_cdna_k15", "-k", "15"]
            subprocess.run(salmon_index_command)

            # Run Salmon on the trimmed FASTQ files
            salmon_command = ["salmon", "quant", "-i", "t_indxs/pao1_cdna_k15", "-l", "ISR", "-r", os.path.join(run_number_folder, trimmed_fastq_files[0]), "--validateMappings", "-o", sra_study_folder]
            print("Running Salmon on the trimmed FASTQ files in run number:", run_number)
            subprocess.run(salmon_command)
            # Move the quant.sf file to Study folder
            os.rename(os.path.join(sra_study_folder, "quant.sf"), os.path.join(sra_final_dir, f"{sra_study}_{trimmed_fastq_files[0][:-3]}_quant.sf"))


        elif len(files_in_run) == 2:
            # Run FastQC on each FASTQ file
            for fastq_file in files_in_run:
                fastqc_command = ["fastqc", os.path.join(fastq_directory, fastq_file)]
                print("Running FastQC on:", fastq_file)
                subprocess.run(fastqc_command)
                html_files = [file for file in os.listdir(fastq_directory) if file.endswith(".html")]
                # Move the HTML file to the sra_study_folder
                html_file_path1 = os.path.join(fastq_directory, html_files[0])
                shutil.move(html_file_path1, os.path.join(sra_final_dir, html_files[0]))

            
            # Run Trim Galore on each FASTQ file
            for fastq_file in files_in_run:
                trim_galore_command = ["trim_galore", os.path.join(fastq_directory, fastq_file), "--output_dir", run_number_folder]
                print("Running Trim Galore on:", fastq_file)
                subprocess.run(trim_galore_command)
                txt_files = [file for file in os.listdir(run_number_folder) if file.endswith(".txt")]
                # Move the txt file to the sra_study_folder
                txt_file_path1 = os.path.join(run_number_folder, txt_files[0])
                shutil.move(txt_file_path1, os.path.join(sra_final_dir, txt_files[0]))
  

            # Run Salmon indexing
            salmon_index_command = ["salmon", "index", "-t", "./t_indxs/pao1_cnda.fa.gz", "-i", "./t_indxs/pao1_cdna_k15", "-k", "15"]
            subprocess.run(salmon_index_command)

             # Identify the trimmed FASTQ file(s)
            trimmed_fastq_files = [file for file in os.listdir(run_number_folder) if file.endswith(".fq")]


            # Run Salmon on the trimmed FASTQ files
            salmon_command = ["salmon", "quant", "-i", "t_indxs/pao1_cdna_k15", "-l", "ISR", "-1", os.path.join(run_number_folder, trimmed_fastq_files[0]), "-2", os.path.join(run_number_folder, trimmed_fastq_files[1]), "--validateMappings", "-o", sra_study_folder]
            print("Running Salmon on 2 trimmed FASTQ files in run number:", run_number)
            subprocess.run(salmon_command)
            # Move the quant.sf file to Study folder
            os.rename(os.path.join(sra_study_folder, "quant.sf"), os.path.join(sra_final_dir, f"{sra_study}_{trimmed_fastq_files[0][:-3]}_{trimmed_fastq_files[1][:-3]}_quant.sf"))


    # Remove the fastq folder
    shutil.rmtree(os.path.join(sra_study_folder, "fastq"))


    # Remove the SRA folder
    shutil.rmtree(os.path.join(sra_study_folder, "SRA"))

    # Remove the logs folder
    shutil.rmtree(os.path.join(sra_study_folder, "logs"))

    # Remove the libParams folder
    shutil.rmtree(os.path.join(sra_study_folder, "libParams"))

    # Remove the aux_info folder
    shutil.rmtree(os.path.join(sra_study_folder, "aux_info"))

    # Remove the lib_format_counts.json file
    os.remove(os.path.join(sra_study_folder, "lib_format_counts.json"))

    # Remove the cmd_info.json file
    os.remove(os.path.join(sra_study_folder, "cmd_info.json"))