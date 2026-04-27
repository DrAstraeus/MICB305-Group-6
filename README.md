# MICB 305: Data Science in Microbiology and Immunology Research
### Vitamin D supplementation is correlated with gut microbial diversity, abundance of specific taxa, and predicted functional potential of metabolic pathways in patients with multiple sclerosis
Vivian Fung, Nelson Kung, Benjamin Lee, Ethan Sun, Ivory Zhang

Dept. of Microbiology and Immunology, University of British Columbia, Vancouver, Canada

## Overview
Multiple Sclerosis (MS) is a chronic central nervous system disorder that frequently results in neurological disability, largely due to an inflammation that damages the myelin sheath. Associations between vitamin D intake and expanded disability status scale (EDSS), a metric quantifying MS disease progression, are inconsistent across clinical trials, suggesting that vitamin D’s impact on the progression of MS is currently not well known. We aim to clarify the association between vitamin D and MS by exploring if and how the gut microbiome may contribute to the relationship between oral vitamin D intake and EDSS.
## Research Question
*Is Vitamin D intake associated with changes in microbiome composition and/or function in multiple sclerosis patients?*
## Dataset
Metadata and 16S rRNA amplicon sequencing datasets of MS patients and household healthy controls were derived from [Zhou et al](https://pmc.ncbi.nlm.nih.gov/articles/PMC10143502/). A diet analysis was conducted to measure nutrient intake, including oral vitamin D supplementation (quantified using international units), by participant questionnaire. In total, the study included 1,152 participants, comprising 576 MS patients and 576 healthy controls, and concluded that there were distinct microbial networks between the two groups.
### Key Variables of Interest
| Variable | Definition |
| ----------------- | ----------- |
| `vitamin D (IU)`, `vitamin.D..IU` | Quantified Vitamin D intake measured in IU (international units) |
| `edss` | Quantified neurological disability in multiple sclerosis patients; uses a twenty-point scale from 0.0 (normal) to 10.0 (death due to MS), with half-point increments |
| `disease` | Patient group (categorical variable); MS, Control |
## Methods
Please view the [shell](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/shell.sh) and [R scripts](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/script.R) generated in our study in the `/scripts` directory. All [input files](https://github.com/DrAstraeus/MICB305-Group-6/tree/main/datasets) can be found in the `/datasets` directory.
- QIIME 2 (Amplicon-2025.4)
  - Manifest Importing
  - DADA2 Denoising
  - Silva (release 138) Taxonomic Classification
  - Alpha Rarefaction
- PICRUSt2 
  - Pathway Abundance
- RStudio (2026.01.0, Build 392)
  - Data Wrangling, Visualization
  - [Alpha Diversity Analysis](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/individual_analyses/ES_Alpha_Diversity.R)
  - [Beta Diversity Analysis](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/individual_analyses/ES_Beta_Diversity.R)
  - [Indicator Taxa Analysis](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/individual_analyses/IZ_NK_IndicTaxa.R)
  - [Differential Abundance Analysis](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/individual_analyses/VF_DiffAbun.R)
  - [PICRUSt2 Functional Analysis](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/individual_analyses/BL_PICRUSt2.R)
  - [Spearman's Correlation Heatmap Analysis](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/scripts/individual_analyses/BL_CorrMat.R)
## Results
Please view the main and supplemental [figures](https://github.com/DrAstraeus/MICB305-Group-6/tree/main/figures) generated in our study in the `/figures` directory.
- Vitamin D supplementation is associated with decreased microbiome complexity and evenness.
  - [Fig. 1](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig1.jpeg)
- Vitamin D supplementation has minimal influence on the gut microbiota shift observed in MS patients.
  - [Fig. 2A](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig2A.png)
  - [Fig. 2B](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig2B.png)
- Selective microbial shifts associated with vitamin D show limited links to MS severity.
  - [Fig. 3A](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig3A.png)
  - [Fig. 3B](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig3B.png)
  - [Fig. 3C](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig3C.png)
- The predicted functional potential of vitamin K2 biosynthesis and purine nucleotide salvage is reduced in MS patients on vitamin D supplementation.
  - [Fig. 4A](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig4A.png)
  - [Fig. 4B](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig4B.png)
  - [Fig. 4C](https://github.com/DrAstraeus/MICB305-Group-6/blob/main/figures/Fig4C.png)
## Acknowledgements
We would like to express our deepest gratitude to [Claire Sie](https://github.com/claire0s), our graduate teaching assistant, for guiding our research process, troubleshooting issues with our code, and supporting us at every step of our project. We would also like to thank [Brian Shao](https://github.com/bjzsh) for providing comments on language and clarity for our project proposal and manuscript. Finally, we thank the instructor for MICB 305: Data Science in Microbiology and Immunology Research, [Dr. Avril Metcalfe-Roach](https://github.com/armetcal), who designed the course modules to help us implement useful data science tools for microbiome research. 
