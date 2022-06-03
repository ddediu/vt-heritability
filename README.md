# Code, data and results acompanying the paper "The heritability of vocal tract structures estimated from structural MRI in a large cohort of Dutch twins"

This repository contains the data, code and full results supporting the paper 
Dediu, D., Jennings, E.M., van ’t Ent, D., Moisik, S.R., Di Pisa, G., Schulze, J., de Geus, E.J.C., den Braber, A., Dolan, C.V. & Boomsma, D.I. (2022) *The heritability of vocal tract structures estimated from structural MRI in a large cohort of Dutch twins*.

The repository is structured as follows:

- the file `vt-heritability.Rproj` is the `RStudio` project
- `data` contains the data, structured as follows:

  - `input` contains the input data ("raw") data, as follows:
  
    - `ReadMe.txt`: a readme document
    - `landmarks.csv`: a comma-separated CSV file contaning the description of the landmarks used to define the phenotypic measures
    - `measures-STUB.csv`: is an empty stub that keeps the structure of the missing "true" comma-separated CSV file `measures.csv` which, unfortunately, cannot be made public as it contains potentially identifiable information about our participants; without the "real" file the analysis cannot be run, but this file can be requested from the [Netherlands Twin Register](https://tweelingenregister.vu.nl/information_for_researchers/working-with-ntr-data)
    - `phenotype_names_decriptions.csv`: a comma-separated CSV file contaning the full name, short name, domain, type, and description of the phenotypic measures

  - `intermediate` contains various intermediate files, as follows:
  
    - `data_for_SEM-STUB.csv`: is an empty stub that keeps the structure of the missing "true" comma-separated CSV file `data_for_SEM.csv` which is dervied from the missing `measures.csv` and which contains potentially identifiable information about our participants; is automatically generated from `measures.csv` (when available) running the pre-processing scripts (see below)
    - `SEM_results*.csv` and `lavaan_results.csv` contain the actual results of fitting the SEM model using `OpenMX` and `lavaan`, respectively, and are described in the corresponding `*_explanations.md` `Markdown` documents

  - `final` contains the final results of the analysis, encapsulated in the comma-separated CSV file `heritability_results.csv` and described in the `heritability_results.md` `Markdown` document

- `figures` contains various figures (autmatically or manually generated) used in the paper and the supplemenary information

- `code` contains the various files needed to perform the analysis and draw the plots:

  - there are 4 pre-processing scripts that need to be run in sequence (but only when the "real" `measures.csv` file is available!):
  
    - `1_data_preparation.R`: `R`script for data cleaning and preparation for downstream analysis
    - `2_OpenMx_model.R` and `2_lavaan_model.R`: `R` scripts that fit the Genetic Covariance Structure Model (GCSM) to the data using either `OpenMX` or `lavaan`
    - `3_analysis.Rmd`: the `Rmarkdown`  script that puts everything together and peforms the actual analyses and plots; it generates the next file:
    - `full_analysis_report.html`: the Full Analysis Report, a self-contained `HTML` document that encapsulates all the analyses and plots.
    - three figures used to genrate the full analysis report

