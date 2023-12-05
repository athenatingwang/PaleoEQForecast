# PaleoEQForecast

The computational code for fitting each renewal process model and carrying out model-averaging for the paper: Wang, T., Griffin, J., Brenna, M., Fletcher, D., Zeng, J., Stirling, M., Dillingham, P. and Kang, J. (2023) Earthquake forecasting from paleoseismic records, under review.

All folders can be downloaded and saved in the same main folder. Inside the folder "DataFinal", the users should create another folder named "chronologies_all_final". For example, if you use "PaleoEQ" as the main folder, then you should have the following folder paths:

/PaleoEQ/Rcode/

/PaleoEQ/Results/

/PaleoEQ/Demo/

/PaleoEQ/ResultsDemo/

/PaleoEQ/DataFinal/chronologies_all_final/

The .csv data files can be downloaded from https://github.com/griffij/QuakeRates/tree/master/chronologies_all_final. The .csv data files should be downloaded and saved in the folder "/DataFinal/chronologies_all_final/".  

## The folder "Demo" contains the same .R files as the folder "Rcode". The only difference is that the files in the "Demo" folder use very small values for the MCMC iterations and thus take much less time to run. Due to this reason, if you run the "Demo" files you may encounter error messages when runnning ModelAve.R and ModelAveRetro.R as the MCMC chains are not converged. The ResultsDemo.zip file contains all the sample results for which ModelAve.R and ModelAveRetro.R in the "Demo" folder can run without errors. 


## Software dependencies 
To run the R functions and JAGS code, one need to use R >= 4.3.0, and the following R packages: "R2jags", "lattice", "actuar", "MCMCvis", "rjags", "wesanderson", "maps", "coda", "bayesplot", "ggplot2", "rstanarm", "grid".

## Run time and memory requirement
### "Demo" folder
To run the R functions and JAGS code in the files bptfit.allfaults.R, gammafit.allfaults.R, lnormfit.allfaults.R, poisfit.allfaults.R,and weibfit.allfaults.R in the "Demo" folder with n.iter.mc = 6000, n.burnin.mc = 1000, n.thin.mc = 100, and K = 50 MC samples, it will take less than 30 minutes on a standard computer (each). The memory requirement for running the ModelAve.R and ModelAveRetro.R files using the output from n.iter.mc = 6000, n.burnin.mc = 1000, n.thin.mc = 100, and K = 50 MC samples is < 2G. Instructions are in each .R file.

### "Rcode" folder
To run the R functions and JAGS code in the files bptfit.allfaults.R, gammafit.allfaults.R, lnormfit.allfaults.R, poisfit.allfaults.R,and weibfit.allfaults.R in ##"Rcode"## folder with n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000, it will take more than four weeks on a standard computer (each). The memory requirement for running the ModelAve.R and ModelAveRetro.R files using the output from n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000 is > 64G. It's faster to run these on supercomputers. Instructions are in each .R file.


## Fitting each renewal process to the Monte Carlo samples of the earthquake ocurrence times from all 93 fault segments
## Fitting the BPT renewal process 
### bptfit.allfaults.R in "Rcode" and "Demo" folders
This file sources bptfit.R and fits a BPT renewal process to the Monte Carlo samples of the earthquake occurrence times from each fault segment. The default is to use a loop on a PC (the lines within "for (i in 1:nfault){...}") to fit the model to all 93 fault segments. This will take a long time to finish if n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000 are used, which are the values we used in the manuscript.

If one wants to test if the code runs properly, then use the following three lines to reduce the computational time, which most likely will not result in convergence of the MCMC chains, but can test that the code is running. These are the current values set in bptfit.allfaults.R in the "Demo" folder.
n.iter.mc = 6000;
n.burnin.mc = 1000;
n.thin.mc = 100;
gelman.cutoff = 5

However, after testing the code, for convergence of MCMC chains, a minimum of the following values are suggested.
n.iter.mc = 55000;
n.burnin.mc = 5000;
n.thin.mc = 10;
gelman.cutoff = 1.2

## Fitting the gamma renewal process 
### gammafit.allfaults.R in "Rcode" and "Demo" folders
This file sources gammafit.R and fits a gamma renewal process to the Monte Carlo samples of the earthquake occurrence times from each fault segment. The default is to use a loop on a PC (the lines within "for (i in 1:nfault){...}") to fit the model to all 93 fault segments. This will take a long time to finish if n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000 are used, which are the values we used in the manuscript.

If one wants to test if the code runs properly, then use the following three lines to reduce the computational time, which most likely will not result in convergence of the MCMC chains, but can test that the code is running. These are the current values set in gammafit.allfaults.R in the "Demo" folder.
n.iter.mc = 6000;
n.burnin.mc = 1000;
n.thin.mc = 100;
gelman.cutoff = 5

However, after testing the code, for convergence of MCMC chains, a minimum of the following values are suggested.
n.iter.mc = 55000;
n.burnin.mc = 5000;
n.thin.mc = 10;
gelman.cutoff = 1.2

## Fitting the lognormal renewal process 
### lnormfit.allfaults.R in "Rcode" and "Demo" folders
This file sources lnormfit.R and fits a lognormal renewal process to the Monte Carlo samples of the earthquake occurrence times from each fault segment. The default is to use a loop on a PC (the lines within "for (i in 1:nfault){...}") to fit the model to all 93 fault segments. This will take a long time to finish if n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000 are used, which are the values we used in the manuscript.

If one wants to test if the code runs properly, then use the following three lines to reduce the computational time, which most likely will not result in convergence of the MCMC chains, but can test that the code is running. These are the current values set in lnormfit.allfaults.R in the "Demo" folder.
n.iter.mc = 6000;
n.burnin.mc = 1000;
n.thin.mc = 100;
gelman.cutoff = 5

However, after testing the code, for convergence of MCMC chains, a minimum of the following values are suggested.
n.iter.mc = 55000;
n.burnin.mc = 5000;
n.thin.mc = 10;
gelman.cutoff = 1.2

## Fitting the Poisson process 
### poisfit.allfaults.R in "Rcode" and "Demo" folders
This file sources poisfit.R and fits a Poisson process to the Monte Carlo samples of the earthquake occurrence times from each fault segment. The default is to use a loop on a PC (the lines within "for (i in 1:nfault){...}") to fit the model to all 93 fault segments. This will take a long time to finish if n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000 are used, which are the values we used in the manuscript.

If one wants to test if the code runs properly, then use the following three lines to reduce the computational time, which most likely will not result in convergence of the MCMC chains, but can test that the code is running. These are the current values set in poisfit.allfaults.R in the "Demo" folder.
n.iter.mc = 6000;
n.burnin.mc = 1000;
n.thin.mc = 100;
gelman.cutoff = 5

However, after testing the code, for convergence of MCMC chains, a minimum of the following values are suggested.
n.iter.mc = 55000;
n.burnin.mc = 5000;
n.thin.mc = 10;
gelman.cutoff = 1.2

## Fitting the Weibull renewal process 
### weibfit.allfaults.R in "Rcode" and "Demo" folders
This file sources weibfit.R and fits a Weibull renewal process to the Monte Carlo samples of the earthquake occurrence times from each fault segment. The default is to use a loop on a PC (the lines within "for (i in 1:nfault){...}") to fit the model to all 93 fault segments. This will take a long time to finish if n.iter.mc = 5010000, n.burnin.mc = 10000, n.thin.mc = 1000 are used, which are the values we used in the manuscript.

If one wants to test if the code runs properly, then use the following three lines to reduce the computational time, which most likely will not result in convergence of the MCMC chains, but can test that the code is running. These are the current values set in weibfit.allfaults.R in the "Demo" folder.
n.iter.mc = 6000;
n.burnin.mc = 1000;
n.thin.mc = 100;
gelman.cutoff = 5

However, after testing the code, for convergence of MCMC chains, a minimum of the following values are suggested.
n.iter.mc = 55000;
n.burnin.mc = 5000;
n.thin.mc = 10;
gelman.cutoff = 1.2

### hpd.interval.R
The code in this file is adapted from the HPDinterval.mcmc() function in the "coda" R package to produce the Highest Posterior Density intervals. I modified the HPDinterval.mcmc() function so that the input argument can be a vector.

## Model-averaged forecasts
### ModelAve.R
This file contains the R code for carrying out model-averaged forecasts. To run this file, one needs all the results files saved in the "Results" folder (or "ResultsDemo" if one runs the Demo files). The files in "Results" are very large (over 30G) so not included in this GitHub repository. If you are interested in getting these result files so that you can run this .R file in the "Rcode" folder, please email me. The current "Results" folder contains all the necessary subfolders and some sample result files. You can run the .R code listed above ModelAve.R to produce the result files as well. After that, you can run ModelAve.R to carry out model-averaged forecast.

## Plotting the figures used in the paper 
### PlotsAndRegression.R
This file provides step-by-step code to produce all the figures appeared in the manuscript and the supplementary file. To run this file, you will need the model-averaged forecasts saved in .image files. You can follow the instructions above to obtain these files or send me an email to obtain a copy.

## Retrospective forecasts for assessment of prediction error
The files in the folder "/Rcode/RetroForecast/" are for carrying out retrospective forecasts. They are very similar to the files bptfit.allfaults.R, gammafit.allfaults.R, lnormfit.allfaults.R, poisfit.allfaults.R, weibfit.allfaults.R and ModelAve.R in the folder "/Rcode/". These files all use the values used in the manuscript so it will take a long while to run. If you are only interested in checking if the code runs, you can use the files in "/Demo/RetroForecast/".



