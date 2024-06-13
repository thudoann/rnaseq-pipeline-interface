import streamlit as st
import pandas as pd
import os
from manual_input import process_manual_input, run_nextflow_pipeline
from csv_input import process_csv_input
from multiqc_viewer import display_study_multiqc
from quant_file_viewer import display_merged_quant_files

st.set_page_config(page_title="RNA-seq Pipeline Interface", layout="wide")

st.title("RNA-seq Pipeline Interface")

# Create tabs for navigation
tabs = st.tabs(["Pipeline Interface", "View MultiQC Report", "View Merged Quantitative Expression"])

# Tab 0: Pipeline Interface
with tabs[0]:
    st.sidebar.title("Pipeline Configuration")

    # Option to choose between Compendium Open Data and Customized Data
    data_source = st.sidebar.radio(
        "Select data source:",
        ('Compendium Open Data', 'Customized Data')
    )

    if data_source == 'Compendium Open Data':
        st.sidebar.header("1. Compendium Open Data")
        input_method = st.sidebar.radio(
            "Select input method:",
            ('Manual Input', 'Upload CSV')
        )

        if input_method == 'Manual Input':
            file_type = st.sidebar.radio(
                "Select your file type:",
                ('Single-end', 'Paired-end')
            )
            SRA_study = st.sidebar.text_input("Enter the Study number:")
            run_numbers = st.sidebar.text_area("Enter the run numbers (separated by spaces):")
            run_numbers_list = run_numbers.split() if SRA_study and run_numbers else []
        else:
            uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=["csv"])
            if uploaded_file:
                df = pd.read_csv(uploaded_file)
                st.sidebar.write("File uploaded successfully:")
                st.sidebar.write(df.head())
                if not all(column in df.columns for column in ['SRA_study', 'Run', 'LibraryLayout']):
                    st.sidebar.error("The CSV file must contain 'SRA_study', 'Run', and 'LibraryLayout' columns.")
                    df = None
            else:
                df = None
    else:
        st.sidebar.header("2. Customized Data")
        input_method = st.sidebar.radio(
            "Select input method:",
            ('Manual Input', 'Upload CSV')
        )

        if input_method == 'Manual Input':
            file_type = st.sidebar.radio(
                "Select your file type:",
                ('Single-end', 'Paired-end')
            )
            SRA_study = st.sidebar.text_input("Enter the Study number:")
            fastq_dir = st.sidebar.text_input("Enter the directory containing FASTQ files:")
        else:
            uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=["csv"])
            if uploaded_file:
                df = pd.read_csv(uploaded_file)
                st.sidebar.write("File uploaded successfully:")
                st.sidebar.write(df.head())
                if not all(column in df.columns for column in ['SRA_study', 'LibraryLayout']):
                    st.sidebar.error("The CSV file must contain 'SRA_study', 'LibraryLayout' columns.")
                    df = None
            else:
                df = None

    kmer_size = st.sidebar.number_input("Enter the k-mer size:", min_value=1, value=15)
    transcriptome = st.sidebar.text_input("Enter the directory containing transcriptome file:")
    cpus = st.sidebar.number_input("Enter the number of CPUs:", min_value=1, max_value=os.cpu_count() - 1)
    email = st.sidebar.text_input("Enter your email address:")

    # Provide a button to run the pipeline
    if st.sidebar.button("Run Pipeline"):
        if data_source == 'Compendium Open Data':
            if input_method == 'Manual Input':
                process_manual_input(file_type, SRA_study, run_numbers_list, kmer_size, cpus, transcriptome, email)
            elif input_method == 'Upload CSV':
                process_csv_input(df, kmer_size, cpus, email, transcriptome)
        elif data_source == 'Customized Data':
            if input_method == 'Manual Input':
                if fastq_dir:
                    run_nextflow_pipeline(file_type, SRA_study, kmer_size, cpus, transcriptome, email, fastq_dir)
                else:
                    st.sidebar.error("Please enter the directory containing FASTQ files.")
            else:
                if df is not None:
                    for index, row in df.iterrows():
                        file_type = row['LibraryLayout'].upper()
                        SRA_study = row['SRA_study']
                        fastq_dir = row['fastq-dir']
                        run_nextflow_pipeline(file_type, SRA_study, kmer_size, cpus, transcriptome, email, fastq_dir)
                else:
                    st.sidebar.error("Please upload a CSV file for Customized Data.")

    st.sidebar.header("Instructions")
    st.sidebar.info("""
    1. **Choose the Data Source**: Select either **Compendium Open Data** or **Customized Data**.
       
    2. **For Compendium Open Data**:
       - **Select Input Method**: Choose either **Manual Input** or **Upload CSV**.
         - **If Manual Input**:
           - Choose the file type (Single or Paired).
           - Enter the study number and run numbers.
         - **If Uploading a CSV**:
           - Ensure the CSV file includes columns for:
             - Study Number
             - SRA Numbers
             - File Type (Single or Paired)
    3. **For Customized Data**:
       - **Select Input Method**: Choose either **Manual Input** or **Upload CSV**.
         - **If Manual Input**:
           - Enter the directory path containing the FASTQ files.
         - **If Uploading a CSV**:
           - Ensure the CSV file includes columns for:
             - Study Number
             - File Type (Single or Paired)
             - Directory path

    4. **Set the k-mer Size**: Enter the k-mer size for indexing (it must be an odd number).

    5. **Set the number of CPUs**: Enter the number of CPUs to be used for the pipeline.

    6. **Enter your email address**: Provide your email address to receive a summary of the pipeline execution.

    7. **Start the Process**: Click on the "Run Pipeline" button to begin.
    """)

# Tab 1: View MultiQC Report
with tabs[1]:
    display_study_multiqc()

# Tab 2: View Merged Quant Files
with tabs[2]:
    display_merged_quant_files()
