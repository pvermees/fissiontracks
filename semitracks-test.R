rm(list=ls())
graphics.off()
source('semitracks.R')

example <- 3

if (example==1){
    fname <- 'track-lengths/fct_c_semi.csv'
    confined <- FALSE
    skip <- 0
} else if (example==2){
    fname <- 'track-lengths/UM10-10 Brogo 3D Semi Tracks with Dper.csv'
    confined <- FALSE
    skip <- 2
} else if (example==3){
    fname <- 'track-lengths/UM10-10 Brogo 3D Confined Lengths.csv'
    confined <- TRUE
    skip <- 2
}

# open the data file
dat <- read.data(fname,confined=confined,cutoff=3,skip=skip)
fit <- invert(dat$l,confined=confined,ncomp=3)
mod <- forward(fit$P,fit$M,fit$S)
plotModel(mod,dat)
