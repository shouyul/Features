## Overview
**Code repository for:**
*A multivariate investigation of visual word, face, and ensemble processing: Perspectives from EEG-based decoding and feature selection*
Psychophysiology, 2020, 57:e13511
[DOI: 10.1111/psyp.13511](https://doi.org/10.1111/psyp.13511)  
Authors: Dan Nemrodov\*, Shouyu Ling\*, Ilya Nudnou, Tyler Roberts, Jonathan S. Cant, Andy C. H. Lee, Adrian Nestor  
(\*Equal contribution)

This repository contains the MATLAB scripts and analysis pipeline used in the study, which applied **EEG-based decoding** and **recursive feature elimination (RFE)** to investigate spatiotemporal and frequency-domain features underlying identity-level visual processing across three categories:

- **Words** (visual word recognition)
- **Faces** (single identity recognition)
- **Face ensembles** (crowds of faces)

Our approach:
- Extracted **time-domain** and **frequency-domain** EEG features
- Ranked feature diagnosticity using **SVM-RFE**
- Assessed **within-** and **cross-participant** stability of features
- Evaluated the boost in classification accuracy and dimensionality reduction achieved through feature selection

---

## Key Findings

1. **Commonalities**  
   - Word and face decoding both relied on **bilateral occipitotemporal channels** and activity peaking around **150–450 ms**.
2. **Differences**  
   - Ensemble decoding engaged **central channels**, with later peaks (~350–450 ms).
   - Alpha-band features were more prominent for **words and ensembles** than for faces.
3. **Feature Selection Benefits**  
   - Boosted classification accuracy, especially for face decoding in the frequency domain.
   - Reduced feature sets by **>80%** in most cases without performance loss.
