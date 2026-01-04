**Fluorescence microscopy colocalization analysis** quantifies the spatial overlap between two or more fluorescent signals in an image in order to assess whether different molecules, structures, or labels are spatially associated.
This is commonly visualized by overlaying fluorescence channels (e.g., red and green, with yellow indicating overlap) and quantified using statistical coefficients such as Pearson’s Correlation Coefficient (PCC) and Manders’ Colocalization Coefficients (MCC).

This repository provides a simple Python pipeline (inspired by FIJI) for performing colocalization analysis on multichannel fluorescence microscopy images.

The pipeline:
- Loads multichannel TIFF images
- Applies thresholding (e.g., Otsu or manual)
- Computes commonly used colocalization metrics
- Generates multiple plots to help interpret and validate the results visually

Metrics being computed:
- [x] **Pearson's Correlation Coefficient (PCC)**
  - Measures linear correlation between intensities of two channels.
    
- [x] **Manders’ Coefficients (M1, M2)**
- Quantify the fraction of signal in one channel that overlaps with the other:

- [x] **Overlap Coefficient**
Measures normalized intensity overlap between channels.

- [x] **Li’s Intensity Correlation Quotient (ICQ)**
- Indicates whether pixel intensities vary synchronously, randomly, or segregated.

- [x] **Cytofluorogram Regression Parameters**
- Computes slope and intercept from the intensity–intensity scatter plot.

Visual outputs:
- [x] Individual red and green channel images
- [x] Merged red and green channels
- [x] Thresholded channel images
- [x] Colocalization mask
- [x] Cytofluorogram plot
- [x] Intensity profiles
- [x] Colocalization intensity maps

Observation:
- This pipeline assumes red and green channels colocalization
- Thresholding strategy can significantly affect results

**TODO / Planned Improvements:**
- [] Costes Thresholding
- [] Background subtraction
- [] Morphological operators
- [] Thresholding Sensitivity Analysis
- [] Metric validation and comparison with existing tools (FIJI/Coloc2/BIOP)


**Contributions are welcome and encouraged.**
This repository is intended to be an open, evolving reference implementation for fluorescence microscopy colocalization analysis. Improvements, bug fixes, extensions, and validation experiments are all appreciated. Please feel free to open an issue to discuss ideas, questions, or potential improvements before submitting a pull request.
  
