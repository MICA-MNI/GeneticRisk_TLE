# Polygenic risk effects of temporal lobe epilepsy on brain morphology in typical development

## Description
**Authors:**   
Alexander Ngo: alexander.ngo@mail.mcgill.ca   
Boris C. Bernhardt: boris.bernhardt@mcgill.ca (corresponding author)   

**Cite:**   
**DOI:**   
**Preprint:** https://www.biorxiv.org/content/10.1101/2025.01.17.633277v1   

This repository contains the code to follow the workflow for our imaging-genetic analysis

## Repository content

 ```
├── README.md
│
├── data
│   ├── raw
│   ├── processed
│   └── results
│       ├── 01_geneticCorrelation
│       ├── 02_epicentreMapping
│       ├── 03_atrophyAssociation
│       ├── 04_epicentreAssociation
│       ├── s01_subcorticalCorrelation
│       ├── s02_thresholdConsistency
│       ├── s03_epicentreConsistency
│       ├── s04_atrophyConsistency
│       ├── s05_diseaseEpicentreConsistency
│       ├── s06_psychiatryAtrophySpecificity
│       └── s07_psychiatryEpicentreSpecificity
│
├── requirements.txt   
│
└── src
    ├── datasets.py
    ├── analyses
    │   ├── 01_geneticCorrelation.py
    │   ├── 02_epicentreMapping.py
    │   ├── 03_atrophyAssociation.py
    │   ├── 04_epicentreAssociation.py
    │   ├── s01_subcorticalCorrelation.py
    │   ├── s02_thresholdConsistency.py
    │   ├── s03_epicentreConsistency.py
    │   ├── s04_atrophyConsistency.py
    │   ├── s05_diseaseEpicentreConsistency.py
    │   ├── s06_psychiatryAtrophySpecificity.py
    │   └── s07_psychiatryEpicentreSpecificity.py
    └── plots
        ├── f01_geneticCorrelation.py
        ├── f02_epicentreMapping.py
        ├── f03_atrophyAssociation.py
        ├── f04_epicentreAssociation.py
        ├── sf01_subcorticalCorrelation.py
        ├── sf02_thresholdConsistency.py
        ├── sf03_epicentreConsistency.py
        ├── sf04_atrophyConsistency.py
        ├── sf05_diseaseEpicentreConsistency.py
        ├── sf06_psychiatryAtrophySpecificity.py
        └── sf07_psychiatryEpicentreSpecificity.py

```

--------
