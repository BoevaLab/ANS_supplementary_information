# EMT Scoring Experiments - README

This folder contains the notebooks required to recreate the EMT (Epithelial-Mesenchymal Transition) signature published in the manuscript.

## Prerequisites

- **Annotations and Preprocessed Data:**  
  Please ensure that all necessary cell annotations and the preprocessed datasets have been downloaded beforehand, c.f., see `README.md` in the root folder.
  You may also need to download the EMT barcodes if required by the analysis.

## Steps to Reproduce the EMT Signatures

1. **Run Cancer-Type Specific Notebooks:**  
   For each cancer type (see subfolders), execute the `get_cancer_emt_cells.ipynb` notebook to identify EMT cells.

2. **Find EMT-Specific Signatures:**
   - For ESCC (Esophageal Squamous Cell Carcinoma):  
     Run `ESCC/find_cancel_emt_signature_ESCC.ipynb` to identify an EMT-specific ESCC signature.
   - For LUAD (Lung Adenocarcinoma):  
     Run `LUAD_xing/find_cancel_emt_signature_LUAD.ipynb` to identify an EMT-specific LUAD signature.

3. **Combine and Refine Signatures:**  
   Run `ESCC/union_ESCC_and_LUAD_specific_EMT_signature_and_refine_on_ESCC.ipynb` to create a union of the ESCC and LUAD EMT signatures and further refine it on ESCC data.

4. **Evaluate the New Signature:**  
   Use `evaluation_LUNG_ESCC_signatures.ipynb` to evaluate the performance of the newly derived signature.

5. **Correlate Signature with Histotypes (e.g., for Figure 3c):**  
   To correlate the EMT signature with histotypes, run the scripts provided in the `correlation_histotypes_R` folder.

---

**Note:**  
All analyses assume that the required data files are present in the expected locations. Please refer to the individual notebooks for additional details and requirements.
