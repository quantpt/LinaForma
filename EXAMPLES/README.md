## ICSV_forward_model.csv
This file contains the forward models for each grid point, which can be created in any forward modelling software (Perple_X, Theriak-Domino, MAGEMin etc.). Make sure the first two columns contain the temperature and pressure points (other any other 2 variables of interest e.g., aH2O-T) that define the grid. The next columns should contain the measured variables for your system of interest (e.g., XAlm). 

## Measurement files
You can either give your measurements as a mean + standard deviation (e.g., InputA.csv) or give all the measurements and the code will calculate the required metrics for you (e.g., InputB.csv). Giving all the measurements also allows you to run non-parametric bootstrap re-sampling.
