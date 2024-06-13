import streamlit as st
import os
from utils import display_multiqc_report

def display_study_multiqc():
    st.title("View MultiQC Report")
    
    # Load available study numbers (assuming each study has a directory in the "Output" folder)
    output_dir = "Output"
    studies = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
    
    if studies:
        selected_study = st.selectbox("Select Study Number", studies, key="multiqc_study_select")
        
        if selected_study:
            multiqc_report_path = os.path.join(output_dir, selected_study, "multiqc_report.html")
            multiqc_report = display_multiqc_report(multiqc_report_path)
            if multiqc_report:
                st.components.v1.html(multiqc_report, height=800, scrolling=True)
            else:
                st.error("MultiQC report not found.")
    else:
        st.error("No studies found.")
