This repository contains all code associated with the manuscript "Directionality of Transcriptional Regulatory Elements" (preprint: https://www.biorxiv.org/content/10.1101/2024.12.01.621925v1).

In this project, we systematically compared the differences between unidirectional and divergent transcriptional regulatory elements (TREs) detected by PRO-cap assay in terms of architecture, function, and evolution. Furthermore, we identified transcription factors with dual roles in transcription initiation, which act in a position-dependent manner and influence the directionality of TREs.

# Installation & Setup
- For most packages, you can create a conda environment using `environment.yml`. For DESeq2-related R packages, use `R.yml` to create a separate conda environment. For tools not listed in `.yml` files, check the source links in the `Tools` section of the `data_prerocessing/2.Other_resources.ipynb` notebook.
- The base directories for input and output data are listed at the top of each notebook. You can modify these hardcoded paths to match your own setup.

# Data Source & Preprocessing
- For PRO-cap data generated in this study, you can download and preprocess the data by following the instructions in notebooks `1-1` to `1-4` (`PROcap_data_preparation` through `TRE_processing`) in the `data_prerocessing` folder.
- For data obtained from public sources, download links and preprocessing steps are provided in the `data_prerocessing/2.Other_resources.ipynb` notebook. 

# Figure Reproduction
- To reproduce the figures in the manuscript, run the notebooks located in the corresponding folders. The code is organized by figure number. For instance, `Fig.1_supps` contains the scripts associated with Fig. 1 and its related figures.

# Contributors
- You Chen (https://github.com/cybluetree/)
- Alden K. Leung (https://github.com/aldenleung/)