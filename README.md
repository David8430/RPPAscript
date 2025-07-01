# RPPAscript

A script in R that uses the core of the RPPASPACE package, but has been modified to accept slide input with sample/dilution replicates.
The script also uses normalisation with general additive models for covariance with spot spacial location and sample protein concentration (determined on a separate slide). Based on the publication of (but may not be exatly the same as) Leanne de Koning et al. 2012.
The script has very little error handling built into it, so reqires proper input.

The app can be run by typing in the R console:
library(shiny)
runApp("filepath\\Shiny_UI.R")

The following packages are required to be installed including their dependencies:
RPPASPACE
mgcv
dplyr
tibble
shiny

The scipts need to be in the proper structure:
AppFolder/Shiny_UI.R
AppFolder/R/customRPPASPACEfunctions.R
AppFolder/R/IterativeProcess.R
AppFolder/R/main.R
AppFolder/R/PostFittingAdj.R
AppFolder/R/PostFittingOutput.R
AppFolder/R/PreFittingDataHandling.R

The slide input information also has to be in te proper structure: (The ProjectFolder is a flexible input, can be called anything else.)
ProjectFolder/in/FCF.csv
ProjectFolder/in/Akt.csv
ProjectFolder/in/FGFR2.csv
etc.
The rest of the folders will be created by the script.

For the directory input the ProjectFolder path needs to be given (example c:\user\paul\Desktop\20210527)

The script can only handle a single series of data, meaing any number of unique antigen/antibody slides and ONE respective FCF protein slide.
Only samples with a matching identifier in the FCF slide AND properly evaluable on both will be returned.

The script does the following:
filters out all spots that are less than the 99 percentile of the blanks (dots without dilution values)
discards all samples that have less than 3 remaining dilutions
optionally flips the source well plate and associated samples the other way around
calculates relative concentration for the slide with cobs fitting before and after spatial correction.
makes intensity corrections based on the spatial covariance of deviation from the response curve
each antigen/antibody slide is corrected for protein covariance
the final protein corrected results are put out as multipliers (transformed back from the logarithmic space)
