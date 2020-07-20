# Source winn.R
source("/an/vital/vital200/WiNN/winn.R")

## The data directory
data.dir <- "/an/vital/vital200/WiNN/data"

## Set study id
# stdy.id <- "synth_linear"
# stdy.id <- "synth_horiz"
stdy.id <- "synth_harm"
# stdy.id <- "VITAL200"

## Load dataset to correct
# met.dat <- get(load(paste(data.dir, "/synthetic/synth_linear.RData",sep="")))
# met.dat <- get(load(paste(synt.dir, "/synthetic/synth_horiz.RData" ,sep="")))
met.dat <- get(load(paste(data.dir, "/synthetic/synth_harm.RData"  ,sep="")))
# met.dat <- get(load(paste(data.dir, "/synthetic/vital_bsln_uncorrected.RData",sep="")))

## Apply correction
correction <- winn(met.dat=met.dat, stdy.id=stdy.id)

