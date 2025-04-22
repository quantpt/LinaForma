 <p align="left">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/d2577b0a12c168a8a8fe5a055eeb452f473757e5/images/logo_black.jpg", width="30%">
</p>

## Overview
These five codes perform a grid-search inversion, accompanied by bootstrap re-sampling, to determine which Pressure-Temperature (P-T) conditions best fit the rock of interest. The bootstrap re-sampling allows the user to determine the uncertainty in this result, as well as the sensitivity of the result to the uncertainty of the input variables. An example of the ouputs can be found in [EXAMPLES](https://github.com/TMackay-Champion/LinaForma/tree/8486dc1820e7d5363f01476148a69ec186ac12be/EXAMPLES).

## Data Inputs
All of the codes require the same two inputs:
1) **forward_model.csv** = this CSV file contains the forward models for selected variables (e.g., XMg in biotite, vol% garnet) created in software like THERIAK-DOMINO. These files can be created using scripts E2 or E3 (see [TheriakDomino_tools](https://github.com/TMackay-Champion/LinaForma/tree/8486dc1820e7d5363f01476148a69ec186ac12be/TheriakDomino_tools)). The codes use the same input format as script E4.
2) **measurements.csv** = this CSV file contains the measured values matching each of the variables in the forward model file. It is important that the measurements are in the order given the the forward model file. There are two accepted formats for this file: A) list all the measurements for each parameter. These may be individual point measurements for example (e.g., [InputA](https://github.com/TMackay-Champion/LinaForma/blob/fc10a0389be120343103fd7d5d064678d722b435/EXAMPLES/InputA.csv)). B) provide the mean and standard deviation of each parameter (e.g., [InputB](https://github.com/TMackay-Champion/LinaForma/blob/fc10a0389be120343103fd7d5d064678d722b435/EXAMPLES/InputB.csv)). 


## L0_isopleths.m
This code allows the user to ascertain which regions in P-T space coincide with the observed values for each parameter of interest. Zones of intersecting isopleths will be plotted.

<details>
<summary> L0 script input </summary>

 % ====== Data ======\
**model = '?'**\
 This is the CSV file for the forward models (e.g., inputs/forward_model.csv).\
**measurements = '?'**\
This is the CSV file for the measurements (e.g., 'inputs/measurement_distributions.csv').

% ====== Data type ======\
**raw = ?**\
This is the type of format you have used for your input measurements file. If you have used InputA (see above), then raw = 1. If you have used InputB, then raw = 0.

% ====== Range of values (only applicable if raw = 0 ======\
**sd = ?**\
This parameter is only applicable if you have used InputA, and controls the range of isopleth values plotted\
(i.e., range = MEAN +/- sd * STANDARDDEVIATION).

% ====== PLOTS ======\
% PLOT 1 = percentage overlap plot\
**all1 = ?**\
% Do you want to plot all of the variables? 1 = YES, 0 = NO.\
**columns1 = [?]**\
If you have selected all1 = 0, which column of the input measurements (i.e., variables) do you want to plot?\

% PLOT 2 = individual isopleths\
**all2 = ?**\
Do you want to plot all of the variables? 1 = YES, 0 = NO.\
**columns2 = [?]**\
If you have selected all2 = 0, which column of the input measurements (i.e., variables) do you want to plot?\
</details>

<details>
<summary> L0 script output </summary>
 
The code outputs three figures: 
1) **Percentage overlap**. This plot shows the regions in P-T space which have the greatest percentage of overlapping variables. 
2) **Isopleths**. This plot shows which regions in P-T space coincide with the observed values for each parameter of interest. Different variables are ascribed different colours.
3) **Overlapping contours**. This plot shows the contours for each parameter and the overlapping areas in P-T space for the measured values.
</details>


## L1_inversion.m
This code performs a bootstrap re-sampling of the observation data, and a grid-search inversion to find the best-fit solution for each bootstrap resample. These are combined to give the user a best-fit pressure-temperature estimate and associated uncertainty for the system of interest. Diagnostic metrics are given with each run so the user can assess the quality of the inversion data-fit.

The grid-search inversion works by calculating the misfit between the observed data and the forward model at each P-T point in the forward model CSV file dataset. The best-fit solution is the P-T point with the lowest misfit value. It is important that the grid spacing is dense enough to ensure the minimum misfit is not missed between the grid spaces. 

<details>
<summary> Diagrammatic representation of the grid-search inversion workflow </summary>
 <p align="center">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/3aaf53b7526049c99e900da48fb3ca8a4db37272/images/L_gridsearch.png", width="90%">
</p>
</p>
</details>

<details>
<summary> Diagrammatic representation of the bootstrap re-sampling workflow </summary>
 <p align="center">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/3aaf53b7526049c99e900da48fb3ca8a4db37272/images/L_bootstrap.png", width="90%">
</p>
</details>

<details>
<summary> L2 script input </summary>
 
% ====== Data ======\
**model = '?'**. As above.\
**measurements = '?'**. As above.

% ====== Data type ======\
**raw = ?**. As above.

% ====== Bootstrapping variables ======\
**bootstrap_type = ?**\
The type of bootstrapping. Parametric = 1, non-parametric = 0. Only parametric is available if raw = 0.\
**it = ?**\
The number of re-samples you want to take.

% ====== PLOTS ======\
**confidence_level = ?**\
The 2D confidence level for the ellipse plotted in Figure 2.\
**boxplots = ?**\
Do you want boxplots or histograms? 1 = boxplot, 0 = histogram.\
**plot_type = ?**\
% What type of plot do you want in Figure 2? 1 = contour plot, 0 = heatmap.\
**T_bins = ?**\
Number of temperature bins in 2D histogram (Figure 2, 4).\
**P_bins = ?**\
Number of pressure bins in 2D histogram (Figure 2, 4).\
</details>

<details>
<summary> L2 script output </summary>
 
The code outputs four figures:

1) a grid showing the extent and resolution of the forward models.
2) the grid-search solution with uncertainty analysis.
3) a plot showing all of the best-fit solutions overlain on the overlapping contour plot of *L0_isopleths.m* script.
4) a 2D histogram of the best-fit solutions.
</details>

## Extra_residuals.m
This script allows the user to examine the difference between the forward model predicitions and the observed data at chosen P-T points.
This can be used to check how well different variables match the best-fit solution. When accompanied by textural evidence or large enough datasets, this process could be used to examine disequilibrium and/or model error.

<details>
<summary> L3 script input </summary>
 
% ====== Data ======\
**model = '?'**. As above.\
**measurements = '?'**. As above.

% ====== Data type ======\
**raw = ?**. As above.

% ====== Select P-T point for forward model data ======\
**T_best = ?**.\
Select the T point of the forward model to which the variables will be compared.\
**P_best = ?**.\
Select the P point of the forward model to which the variables will be compared. Units = bars.
</details>

<details>
<summary> L3 script ouput </summary>
This code outputs boxplots for each variable showing the distribution of observations and the forward model predicted value.
</details>


## Extra_sensitivity.m
This script examines how sensitive the best-fit solutions are to uncertainty in the observations. To do this, the code bootstrap re-samples one variable at a time while keeping all the other variables constant. The resulting variation in the best-fit solutions can then be directly attritubted to the uncertainty on that particular variable. 

<details>
<summary> Diagrammatic representation of the sensitivity workflow </summary>
 <p align="center">
<img src="https://github.com/TMackay-Champion/LinaForma/blob/3aaf53b7526049c99e900da48fb3ca8a4db37272/images/L_sensitivity.png", width="90%">
</p>
</details>

<details>
<summary> L4 script input </summary>
 
% ====== Data ======\
**model = '?'**. As above.\
**measurements = '?'**. As above.

% ====== Data type ======\
**raw = ?**. As above.

% ====== Bootstrapping variables ======\
**bootstrap_type = ?**. As above.\
**it = ?**. As above.\

% ====== Select P-T point for sensitivity analysis ======\
**T_best = ?**\
Select the best-fit T point at which the sensitivity will be computed.
**P_best = ?**\
Select the best-fit P point at which the sensitivity will be computed.
</details>

<details>
<summary> L4 script output </summary>
This code outputs two "tornado" plots, one for temperature and one for pressure.
These plots display how the variation in a particular variable influences the best-fit solutions relative to a given best-fit solution (ideally the output of the *L2_inversion.m* script).
</details>


## FAQs
<details>
<summary> How many input variables should I use? </summary>
A grid-search is a non-linear inversion. The problem is overdetermined because the number of observations is in excess of of the number of model parameters (T, P). Each data point provides a constraint on the possible solution. By incorporating multiple constraints, overdetermined problems can identify and compensate for errors in the different variables.This often results in a higher precision estimate than could be achieved by the individual variables alone. 

This is important when considering correlated variables. In an error-free system, it would be reasonable to remove all correlated variables. However, each variable in a petrological system is associated with a different level of error and the primary cause of this error will vary between different variables. One can readily imagine a situation in which two highly correlated variables result in different P-T estimates due to petrological or model error. Considering the two variables together allows the user to determine the most appropriate solution.
As such, we deem it acceptable to use variables which are predicted to be highly correlated in the forward model. However, phases which are common to one mineral should only be used together if a degree of freedom remains. For example, e.g., the composition of plagioclase can be described by XAn and XAb, or garnet by XAlm, XGrs, XPrp, and XSpss. To maintain a degree of freedom you should use one fewer than the total number of parameters defining the mineral system.  
</details>

<details>
 <summary> What are the different types of bootstrap re-sampling? </summary>
Bootstrap resampling refers to random re-sampling of the original dataset, with replacement. This re-sampling is used to examine the uncertainty of the best-fit solution. Uncertainty may be introduced through model error, mesasurement error, disequilibrium and other such processes. Using bootstrapping, we can examine the effects of these errors on our P-T solution. The codes allow two different styles of re-sampling for each variable.
 
 **(1) Parametric bootstrapping (option 1)**: this assumes that the data follow a specific parametric distribution, such as the normal distribution. 
Bootstrap samples are chosen by randomly drawing observations from the 
distribution with replacement. LinaForma assumes a normal distribtuion. If the measurement data is in input style A, the distrubution is calculated from all the variable measurements. If input style B is used, the distrubution is calculated using the given mean and standard deviation.
 
 **(2) Non-parametric bootstrapping (option 0)**: this option re-samples the original data multiple times with replacement and calculates the mean from each set of samples. These mean values are used in the grid-search. This approach is only possible if the input style A is used.
</details>
 
<details>
<summary> What bootstrapping method should I use? </summary>
 
 The method of bootstrapping depends on your assumptions surrounding the sources of error in the system. Instead of making explicit assumptions about the model, non-parametric bootstrapping focuses solely on the observed data and its properties. This error can be examinied for each variable using the *L1_error.m* script. As such, we deem this bootstrap method to be most appropriate if we assume that the primary source of error is analytical and/or related to disequilbrium, geological uncertainty etc. In this case, the error associated with the observations is greater than associated model error. 

However, in some cases the primary source of error may be model error. In this case, the parametric bootstrap option may be most suitable as it allows the user to select an appropriate mean and standard deviation. A standard deviation should be chosen which allows for a suitable degree of temperature and pressure uncertainty for the particular variable of interest. This can be chosen with the help of the *EXTRA_synthetic_variation.m* script.
</details>










