# Polygenic risk effects of temporal lobe epilepsy on brain morphology in typical development

## Description
**Authors:**   
Alexander Ngo: alexander.ngo@mail.mcgill.ca   
Boris C. Bernhardt: boris.bernhardt@mcgill.ca (corresponding author)   

**Cite:**   
**DOI:**   
**Preprint:**   

This repository contains the code to follow the workflow for our imaging-genetic analysis

## Repository content
 
 ```
├── README.md
│ 
├── data                <- Datasets required for analysis
│   ├── raw
│   ├── processed
│   └── results
│
├── figures             <- Generated necessary plots and graphics for each figure
│   ├── figure1
│   ├── figure2
│   ├── figure3
│   ├── supplementary1
│   └── supplementary2
│
├── requirements.txt   
│
└── src                 <- Source code
    ├── datasets.py      <- Code to load data
    ├── analyses         <- Code to run data analysis
    │   ├── 01_geneticCorrelation.py
    │   ├── 02_epicentreMapping.py
    │   ├── 03_casecontrolAssociation.py
    │   ├── s01_thresholdConsistency.py
    │   └── s02_casecontrolConsistency.py
    └── plots            <- Code to create visualizations
        ├── f01_geneticCorrelation.py
        ├── f02_epicentreMapping.py
        ├── f03_casecontrolAssociation.py
        ├── sf01_thresholdConsistency.py
        └── sf02_casecontrolConsistency.py

```

--------
