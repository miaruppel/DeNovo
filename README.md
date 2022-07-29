# De Novo Peptide Identification - MassSpec
### Important Scripts:
1. DeNovo_main 
    - Need 3 files: mzML, log, realtime search
    - Analysis of MS2
      - How many fragments are missing and at which locations?
    - Integrated R script (fragmentStatistics.R) for data manipulation
    - Matches up MS2 and MS3 scans 
    - Finds which fragments were found in MS3

2. analyzeMS3Data
    - Imports DeNovo_main 
    - Closer analysis of MS3 
      - Check false positive rate 
      - Finds which fragments (that were misisng in MS2) were found in MS3
      - Bar plots for data visualization 
  
3. Pretrained_prosit
    - Data organization intended for pretrained prosit model (HLA-CID)


### Work in Progress:
1. DeNovo_formatData
    - Format data for training (one experiment file)

2. DeNovo_formatData_multiple

4. All prosit helper scripts 
    - match.py, model.py, tensorize.py, etc
