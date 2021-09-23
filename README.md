# bevacizumab-response
A complete directory and sample heirarchy of reproducable scripts for analysis of PDX glioblastoma samples to elucidate genetic contributors to differential response to bevacizumab therapy.

## Metadata
doi: <br>
Author(s): Roshan Lodha<sup>a,b</sup>, Candece Gladson<sup>b</sup> <br>
<sup>b</sup>Lerner College of Medicine, Cleveland Clinic, 9501 Euclid Avenue, Cleveland, Ohio 44195
<sup>b</sup>Lerner Research Institute, Cleveland Clinic, 9500 Euclid Avenue, Cleveland, Ohio 44195


## Analysis of Glioblastoma PDX samples 
The following heirarchy allows for one-shot analysis starting with post-read mapped data. 
```bash
├── raw
│   ├── blood_vessel_geneset.txt
│   ├── counts.csv
│   ├── counts.tabular
│   ├── studydesign.txt
│   ├── prefiltergbmexpr.csv
│   ├── prefilterstudydesign.txt
├── tables
│   ├── tables.R
├── deseq
│   ├── 
├── gsea
│   ├── 
├── dim
│   ├── 
├── cBioPortal
│   ├── data
│   │   ├── CHRNA7_survival_cbioportal.txt
│   ├── plots
│   │   ├── 
│   ├── cBioPortal.R
├── bev_analysis.R
├── README.md
└── .gitignore
```
Specific parameters and version control are detailed in the methods section the paper. 
