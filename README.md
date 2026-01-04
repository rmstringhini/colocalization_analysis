Fluorescence microscopy colocalization analysis quantifies the spatial overlap between two or more fluorescent signals in an image in order to assess whether different molecules, structures, or labels are spatially associated.
This is commonly visualized by overlaying fluorescence channels (e.g., red and green, with yellow indicating overlap) and quantified using statistical coefficients such as Pearson’s Correlation Coefficient (PCC) and Manders’ Colocalization Coefficients (MCC).

This repository provides a simple Python pipeline (inspired by FIJI) for performing colocalization analysis on multichannel fluorescence microscopy images.

The pipeline:
- Loads multichannel TIFF images
- Applies thresholding (e.g., Otsu or manual)
- Computes commonly used colocalization metrics
- Generates multiple plots to help interpret and validate the results visually

Metrics being computed:
**Pearson's Correlation Coefficient (PCC)**
  - Measures linear correlation between intensities of two channels.
    
**Manders’ Coefficients (M1, M2)**
- Quantify the fraction of signal in one channel that overlaps with the other:

**Overlap Coefficient**
Measures normalized intensity overlap between channels.

**Li’s Intensity Correlation Quotient (ICQ)**
- Indicates whether pixel intensities vary synchronously, randomly, or segregated.

**Cytofluorogram Regression Parameters**
- Computes slope and intercept from the intensity–intensity scatter plot.

Visual outputs:
- Individual red and green channel images
- Merged red and green channels
- Thresholded channel images
- Colocalization mask
- Cytofluorogram plot
- Intensity profiles
- Colocalization intensity maps

Observation:
- This pipeline assumes red and green channels colocalization
- Thresholding strategy can significantly affect results

**TODO / Planned Improvements:**
- Costes Thresholding
- Background subtraction
- Morphological operators
- Thresholding Sensitivity Analysis


  
