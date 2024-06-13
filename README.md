# RNA-seq Pipeline Interface

This project provides an RNA-seq pipeline interface using Nextflow for processing RNA-seq data. The interface allows users to configure and run the pipeline through a Streamlit application. The pipeline includes downloading SRA data, extracting FASTQ files, trimming reads, performing quantification, and generating quality control reports.

## Table of Contents
1. [Requirements](#requirements)
2. [Setup](#setup)
3. [Usage](#usage)

## Requirements

- Python 3.6+
- Streamlit
- Nextflow
- Salmon
- Trim Galore
- FastQC
- MultiQC
- SRA Toolkit (`prefetch` and `parallel-fastq-dump`)

Ensure the following tools are installed and accessible in your system's PATH:
- prefetch
- parallel-fastq-dump
- nextflow
- salmon
- trim_galore
- fastqc
- multiqc

## Setup

1. **Clone the repository:**
   ```sh
   git clone https://github.com/thudoann/rnaseq-pipeline-interface.git
   cd rnaseq-pipeline-interface
   ```

2. **Install Python dependencies:**
   ```sh
   pip install -r requirements.txt
   ```

3. **Ensure Nextflow and other tools are installed:**
   Follow the respective installation instructions for Nextflow, Salmon, Trim Galore, FastQC, and MultiQC.

## Usage

1. **Run the Streamlit application:**
   ```sh
   streamlit run app.py
   ```

2. **Configure the pipeline via the Streamlit interface:**
   - Choose between manual input or CSV upload.
   - For manual input, specify the file type (single-end or paired-end), study number, and run numbers.
   - For CSV upload, ensure the CSV contains columns for `SRA_study`, `Run`, and `LibraryLayout`.

3. **Run the pipeline:**
   - Click on the "Run Pipeline" button in the Streamlit interface.
   - The pipeline will download SRA data, process FASTQ files, and execute the RNA-seq analysis.
   - Upon completion, view the MultiQC report directly within the Streamlit app.

