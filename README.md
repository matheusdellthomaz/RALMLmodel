# RALMLmodel
Here you can find the R scripts, RDS files, csv datasets, used in raltegravir (RAL) exposure prediction using machine learning (ML).

Test_combination is the script used for testing different ML algorithms, combination sizes, concentrations at different times, inclusion of covariates and derived predictors from concentrations, developed for predicting RAL AUC from 3 samples.
In this script the algorithms available are xgboost, glm, rf, and svm.
Mapbayaugcomb is the script used for the MAP-BE methodology.

RALscript contains the POPPK models implementation, PK profiles and patients simulation, datasets management, plots, metrics and initial feature selection.

The RDS file named "model_xgboost2021FIM.rds" is the RDS file of the xgboost model developed in this study and the "explainer_externalxgb2021FIM.rds" was used for the shiny webapplication development. 
The RDS can be loaded directly to get predictions from other datasets run the script and get the RAL AUC predictions.

The CSV file named "ralsim2021.csv" is the simulated patients dataset used for training the ML models.
The CSV files named "ralsim2012.csv" is the dataset, for the simulated patients, used for the external validation of the ML models.

The CSV files named "models12hrspairs.csv" and "models12hrstriplets.csv" are the datasets containing the metrics for all the ML algorithms and all sample time combinations of pair and triplets, repectivelly.
