
# Preface -----------------------------------------------------------------

library(TMB)

rm(list=ls())


write('Sys.setenv(PATH=paste(R.home("bin"),Sys.getenv("PATH"),sep=";"))',file="~/.Rprofile",append = TRUE)

## Optionally:
## precompile()
runExample(all=TRUE)


download.file("https://github.com/kaskr/adcomp/archive/master.zip", "adcomp-master.zip")
unzip("adcomp-master.zip")

setwd("adcomp-master/tmb_examples")


# HMM ---------------------------------------------------------------------

source("hmm.R")


# Laplace -----------------------------------------------------------------


source("laplace.R")


# Linear regression -------------------------------------------------------


source("sam.R")


# SDE linear -------------------------------------------------------


source("sde_linear.R")


