import streamlit as st
import pandas as pd
import os

def display_merged_quant_files():
    st.title("Quantitative Expression")
    
    # Load available study numbers (assuming each study has a directory in the "Output" folder)
    output_dir = "Output"
    studies = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
    
    if studies:
        selected_study = st.selectbox("Select Study Number", studies, key="merged_study_select")
        
        if selected_study:
            merged_quant_file_path = os.path.join(output_dir, selected_study, "merged_quant.sf")
            if os.path.exists(merged_quant_file_path):
                st.write(f"Found merged quant.sf file at {merged_quant_file_path}")
                df = pd.read_csv(merged_quant_file_path, sep='\t')
                st.dataframe(df)
            else:
                st.error(f"Merged quant.sf file not found at {merged_quant_file_path}")
    else:
        st.error("No studies found.")
