# PaleoEQForecast

The computational code for fitting each renewal process model and carrying out model-averaging for the paper: Wang, T., Griffin, J., Brenna, M., Fletcher, D., Zeng, J., Stirling, M., Dillingham, P. and Kang, J. (2023) Earthquake forecasting from paleoseismic records, under review.

All folders can be downloaded and saved in the same main folder. Inside the folder "DataFinal", the users should create another folder named "chronologies_all_final". For example, if you use "PaleoEQ" as the main folder, then you should have the following folder paths:

/PaleoEQ/Rcode/

/PaleoEQ/Results/

/PaleoEQ/DataFinal/chronologies_all_final/

The .csv data files can be downloaded from https://github.com/griffij/QuakeRates. The .csv data files should be downloaded and saved in the folder "/DataFinal/chronologies_all_final/".  

## Fitting each renewal process to the Monte Carlo samples of the earthquake ocurrence times from all 93 fault segments
## Fitting the BPT renewal process 
### bptfit.allfaults.R
This file sources bptfit.R and fits a BPT renewal process to the earthquake occurrence times from each fault segment. The default is to use a loop on a PC (the lines within "for (i in 1:nfault){...}") to fit the model to all 93 fault segments. This will take a long time to finish if n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc=1000 are used, which are the values we used in the manuscript.

If one wants to test if the code runs properly, then use the following three lines to reduce the computational time, which most likely will not result in convergence of the MCMC chains, but can test that the code is running. These are the current values set in the bptfit.allfaults.R file.
n.iter.mc = 6000
n.burnin.mc = 1000
n.thin.mc=1
gelman.cutoff <- 2

However, after testing the code, for convergence of MCMC chains, a minimum of the following values are suggested.
n.iter.mc = 55000
n.burnin.mc = 5000
n.thin.mc=10
gelman.cutoff <- 1.2

### gammafit.allfaults.R


### GVPmod0WeibsepZuse.jags 
This file contains the JAGS "model" for the Weibull model in equation (1) in the paper Wang et al. (2022).

### GVPmod1WeibsepZuse.jags
This file contains the JAGS "model" for model M1 in equation (2) in the paper Wang et al. (2022).

### GVPmod2WeibsepZuse.jags
This file contains the JAGS "model" for model M2 in equation (3) in the paper Wang et al. (2022).

### GVPmod3WeibsepZuse.jags
This file contains the JAGS "model" for model M3 in equation (4) in the paper Wang et al. (2022).

### ReadGVPdata.R
This file reads all the data files for the empirical analogues saved in the folder "HoloceneRecords". These data files are obtained from the Smithsonian Institution's GVP catalogue.

##
## When fitting models, make sure to save all .R and .jags files and the two data folders "HoloceneRecords" and "QuaternaryRecords" in the same main folder.

## Fitting models M1, M2, and M3 to the empirical analogues 
### GVPEmpAnalog.R
Run the file GVPEmpAnalog.R, and the MCMC samples of the posterior distributions will be saved in .image files.

## Fitting models M1, M2, and M3 to the statistical analogues with Holocene records, and carry out residual analysis 
### GVPStatsAnalog.R
Run the file GVPStatsAnalog.R, and the MCMC samples of the posterior distributions will be saved in .image files. The residual analysis for each model will be saved as .eps file. If a different format is preferred, such as a .pdf file, just change the "postscript" command to "pdf" and change ".eps" to ".pdf" in the R code.

## Fitting models M1, M2, and M3 to the statistical analogues with Quaternary records (VEI 4+), and carry out residual analysis 
### QuatStatsAnalogVEI4+.R
Run the file QuatStatsAnalogVEI4+.R, and the MCMC samples of the posterior distributions will be saved in .image files. The residual analysis for each model will be saved as .eps file. If a different format is preferred, such as a .pdf file, just change the "postscript" command to "pdf" and change ".eps" to ".pdf" in the R code.

## Fitting models M1, M2, and M3 to the statistical analogues with Quaternary records (VEI 3+), and carry out residual analysis 
### QuatStatsAnalogVEI3+.R
Run the file QuatStatsAnalogVEI3+.R, and the MCMC samples of the posterior distributions will be saved in .image files. The residual analysis for each model will be saved as .eps file. If a different format is preferred, such as a .pdf file, just change the "postscript" command to "pdf" and change ".eps" to ".pdf" in the R code.

## Fitting model M0, a simple Weibull renewal process to the Tongariro GVP Holocene record only 
### GVP-TgrOnly-VEI3+.R
Run the file GVP-TgrOnly-VEI3+.R, and the MCMC samples of the posterior distributions will be saved in .image files. It also provides forecasts from using Holocene record only from Tongariro volcano.

## Fitting model M0, a simple Weibull renewal process to the Tongariro Quaternary record only 
### Quat-TgrOnly-VEI3+.R 
Run the file Quat-TgrOnly-VEI3+.R, and the MCMC samples of the posterior distributions will be saved in .image files. It also provides forecasts from using Quaternary record only from Tongariro volcano.


## Plotting the VEI, posterior distributions, and forecasts used in the paper 
### GRLplots.R
All the results files used in this code are saved in the Results folder.





