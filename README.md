# ILD-BAL-scRNA-Manuscript-2025 Repository

This repository contains the code and data used to generate the figures in the manuscript:

> L.Y. Zhang\*, P.C. Allen\*, _et al._ Transcriptomics of Alveolar Immune Cells Reveals Insight into Mechanisms of Human Pulmonary Fibrosis. _Pending Review._ 

---

## <img src="https://github.com/user-attachments/assets/0ff42912-3def-49d1-8afc-3dc7bb72c401" width="20" height="20"> Table of Contents
- [Overview](#-overview)
- [Installation](#-installation)
- [Usage](#usage)
- [Figure Reproduction Guide](#-figure-reproduction-guide)
- [Data](#-data)
- [Citation](#-citation)
- [Contact](#-contact)

---

## <img src="https://github.com/user-attachments/assets/f451df92-7cda-41eb-906d-85f18984cb96" width="20" height="20"> Overview

This repository provides the source code and data to reproduce the plots and results presented in the manuscript. Each figure in the paper corresponds to a notebook in the `src/` directory.

---

## <img src="https://github.com/user-attachments/assets/8236200d-8a12-45af-997b-890bd235c070" width="20" height="20"> Installation

Clone the repository:

```bash
git clone https://github.com/powellgenomicslab/ILD-BAL-scRNA-Manuscript-2025.git
cd ILD-BAL-scRNA-Manuscript-2025
```

Create a [conda](https://docs.conda.io) environment and install dependencies:

```bash
conda env create -f ild-bal-scrna.yaml
conda activate ild-bal-scrna
```

---

## Usage

Each figure can be reproduced using the associated script or notebook. For example:

For figure 1:

```bash
jupyter notebook src/figure-1.ipynb
```

---

## <img src="https://github.com/user-attachments/assets/a53e141d-5f80-4d1d-9548-bec3d2a148bb" width="20" height="20"> Figure Reproduction Guide

| Figure | Notebook | Notes |
|--------|-------------------|-------|
| Fig. 1 | `src/figure-1.ipynb` | Generates a UMAP plot of all BAL cells post-QC/integration |
| Fig. 2 | `src/figure-2.ipynb` | Subsets the MLMs; generates: 1) Heatmap of Top DEGs between IPF/non-IPF; 2) UMAP of MLMs; 3) Dotplot of DEGs per MLM |
| Fig. 3 | `src/figure-3.ipynb` | MLM Characterization: 1) Pathway; 2) Trajectory; 3) Proportion; 4) DEA of pro-fibrotic genes |
| Fig. 4 | `src/figure-4.ipynb` | Short telomere comparative analysis: 1) Proportion; 2) DEG Volcano |
| Fig. 5 | `src/figure-5.ipynb` | IPF GWAS Integration via scDRS |

---

## <img src="https://github.com/user-attachments/assets/18888156-054a-47ba-9f08-a1c6283912ca" width="20" height="20"> Data

Some figures rely on preprocessed datasets stored in the `data/` directory. Due to size or licensing constraints, some datasets may need to be downloaded manually. See `data/README.md` for instructions.

---

## <img src="https://github.com/user-attachments/assets/f33c77bc-8b93-4d14-8531-ae6e18c84027" width="20" height="20"> Citation

If you use this code or data in your research, please cite:

```
@article{<zhangallen2025>,
  title={Transcriptomics of Alveolar Immune Cells Reveals Insight into Mechanisms of Human Pulmonary Fibrosis},
  author={Zhang, Lai-Ying and Allen, Peter C.},
  journal={pending},
  year={2025},
  doi={pending}
}
```

---

## <img src="https://github.com/user-attachments/assets/9e5b5c7a-bac5-44cb-a568-6e8e03c88017" width="20" height="20"> Contact

For any questions, please contact either:

| Author | Contact |
| :----- | :------ |
| Lai-Ying Zhang   | <a href="mailto:Lai-Ying.Zhang@health.qld.gov.au"><img src="https://github.com/user-attachments/assets/d100bf6e-37cd-458a-aacc-737b26d18b55" width="20" height="20" /></a><a href="https://github.com/laiyingz"><img src="https://github.com/user-attachments/assets/35eceb13-d725-4108-9791-57ca28776e5c" width="20" height="20" /></a> |
| Peter C. Allen | <a href="mailto:p.allen@garvan.org.au"><img src="https://github.com/user-attachments/assets/d100bf6e-37cd-458a-aacc-737b26d18b55" width="20" height="20" /></a><a href="https://petercallen.github.io"><img src="https://github.com/user-attachments/assets/35eceb13-d725-4108-9791-57ca28776e5c" width="20" height="20" /></a> |






 

Otherwise, [raise an issue](https://github.com/PeterCAllen/ILD-BAL-scRNA-Manuscript-2025/issues/new) in the repository and we will address it.
