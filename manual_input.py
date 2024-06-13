import streamlit as st
import subprocess
import os
import time
from utils import run_command, display_multiqc_report, extract_folder_names

def process_manual_input(file_type, SRA_study, run_numbers_list, kmer_size, cpus, email, transcriptome):
    if not SRA_study or not run_numbers_list:
        st.sidebar.error("Please enter both the study number and run numbers.")
    elif kmer_size % 2 == 0:
        st.sidebar.error("The k-mer size must be an odd number.")
    else:
        st.write("## Pipeline Execution")
        with st.spinner("Running script for download FASTQ files..."):
            run_numbers_str = ' '.join(run_numbers_list)
            command = f"python fastqdl.py manual {file_type} {SRA_study} {run_numbers_str}"
            stdout, stderr, returncode, execution_time = run_command(command)
            if returncode != 0:
                st.error("Error running Python script.")
                st.write("### Download FASTQ files Script Error Details")
                st.code(stderr)
            else:
                st.success(f"Download FASTQ files script completed successfully in {execution_time:.2f} seconds.")
                run_nextflow_pipeline(file_type, SRA_study, kmer_size, cpus, email, transcriptome)



def run_nextflow_pipeline(file_type, SRA_study, kmer_size, cpus, transcriptome, email, fastq_dir=None):
    nf_script = f"{file_type.lower().replace('-end', '')}_rnaseq.nf"
    output_dir = os.path.join("Output", SRA_study)
    os.makedirs(output_dir, exist_ok=True)
    with st.spinner(f"An email will be sent when the script finishes."):
        command = f"nextflow run {nf_script} --kmer_size {kmer_size} --study_id {SRA_study} --output_dir {output_dir} --cpus {cpus} --transcriptome {transcriptome}"
        if fastq_dir:
            command += f" --fastq_dir {fastq_dir}"
        #st.write(f"Executing: {command}")  # Print the exact command being run for debugging
        stdout, stderr, returncode, execution_time = run_command(command)
        if returncode != 0:
            st.error("Error running script.")
            st.write("### Script Error Details")
            st.code(stderr)
        else:
            st.success(f"Rnaseq pipeline script ran successfully in {execution_time:.2f} seconds!")
            st.write("## MultiQC Report")
            multiqc_report_path = os.path.join(output_dir, "multiqc_report.html")
            multiqc_report = display_multiqc_report(multiqc_report_path)
            if multiqc_report:
                st.components.v1.html(multiqc_report, height=800, scrolling=True)
            else:
                st.error("MultiQC report not found.")
            
            # Send the summary email
            summary_file = os.path.join(output_dir, "pipeline_summary.txt")
            email_command = f"python send_email.py {summary_file} {email}"
            st.write(f"Executing: {email_command}")
            email_stdout, email_stderr, email_returncode, email_execution_time = run_command(email_command)
            if email_returncode != 0:
                st.error("Error sending summary email.")
                st.write("### Email Script Error Details")
                st.code(email_stderr)
            else:
                st.success("Summary email sent successfully.")

            # Run the merge_quant_files.py script
            try:
                run_numbers_str = extract_folder_names(output_dir)
                #st.write(f"run_numbers_str: {run_numbers_str}")
                merge_command = f"python merge_quant_files.py {SRA_study} {' '.join(run_numbers_str)}"
                #st.write(f"Executing: {merge_command}")
                merge_stdout, merge_stderr, merge_returncode, merge_execution_time = run_command(merge_command)
                if merge_returncode != 0:
                    st.error("Error merging quant.sf files.")
                    st.write("### Merge Script Error Details")
                    st.code(merge_stderr)
                else:
                    st.success(f"Merged quant.sf file created successfully in {merge_execution_time:.2f} seconds.")
            except Exception as e:
                st.error(f"Error running merge_quant_files.py: {e}")
