# RPPAscript

A script in R that uses the core of the RPPASPACE package, but has been modified to accept slide input with sample/dilution replicates.
The script also uses normalisation with general additive models for covariance with spot spacial location and sample protein concentration (determined on a separate slide). Based on the publication of (but may not be exatly the same as) Sylvie Troncale et al. 2012. <br/>
The script has very little error handling built into it, so reqires proper input.

The app can be run by typing in the R console: <br/>
library(shiny) <br/>
runApp("filepath\\Shiny_UI.R") <br/>

The following packages are required to be installed including their dependencies: <br/>
RPPASPACE <br/>
mgcv <br/>
dplyr <br/>
tibble <br/>
shiny <br/>

The scipts need to be in the proper structure: <br/>
AppFolder/Shiny_UI.R <br/>
AppFolder/R/customRPPASPACEfunctions.R <br/>
AppFolder/R/IterativeProcess.R <br/>
AppFolder/R/main.R <br/>
AppFolder/R/PostFittingAdj.R <br/>
AppFolder/R/PostFittingOutput.R <br/>
AppFolder/R/PreFittingDataHandling.R <br/>

The slide input information also has to be in te proper structure: (The ProjectFolder is a flexible input, can be called anything else.) <br/>
ProjectFolder/in/FCF.csv <br/>
ProjectFolder/in/Akt.csv <br/>
ProjectFolder/in/FGFR2.csv <br/>
etc. <br/>
The rest of the folders will be created by the script.

For the directory input the ProjectFolder path needs to be given (example c:\user\paul\Desktop\20210527)

The script can only handle a single series of data, meaing any number of unique antigen/antibody slides and ONE respective FCF protein slide.
Only samples with a matching identifier in the FCF slide AND properly evaluable on both will be returned.

The script does the following: <br/>
filters out all spots that are less than the 99 percentile of the blanks (dots without dilution values) <br/>
discards all samples that have less than 3 remaining dilutions <br/>
optionally flips the source well plate and associated samples the other way around <br/>
calculates relative concentration for the slide with cobs fitting before and after spatial correction. <br/>
makes intensity corrections based on the spatial covariance of deviation from the response curve <br/>
each antigen/antibody slide is corrected for protein covariance <br/>
the final protein corrected results are put out as multipliers (transformed back from the logarithmic space) <br/>

Progress can be tracked on the R console output. When done the windows can simply be closed.
