**"Predicting All-Cause Mortality Using Genetic Data and Biological Aging Markers in the TwinGene Cohort"** 

## Overview 

This project evaluates the predictive performance of biological aging (BA) measures and polygenic risk scores (PRS-mvAge) for all-cause mortality using the Swedish TwinGene cohort. The analysis includes:

- Calculation of BA measures (PhenoAge, KDM, HD, FI, Telomere)
- Construction of PRS for aging-related traits (“mvAge”) using GWAS summary statistics (https://www.nature.com/articles/s43587-023-00455-5)  
- Univariate and multivariable Cox-PH regression models
- Ensemble machine learning models (stacking)
- Subgroup and sensitivity analyses

## Directory Structure

- `TwinGene/scripts/` - R scripts for data analysis and modeling
- `TwinGene/data/` - Input data files (not included in the repository)
- `TwinGene/results/` - Output results such as tables and model outputs
- `TwinGene/figures/` - Generated plots and visualizations

## Requirements

See `requirements.txt` for the list of R packages used in the analysis.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
For questions or contributions, please contact:
•	Dr. Shayan Mostafaei (shayan.mostafaei@ki.se) 
•	Dr. Sara Hägg (sara.hagg@ki.se)

