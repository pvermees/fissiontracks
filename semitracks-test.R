rm(list=ls())
graphics.off()
setwd('~/Documents/Programming/R/fissiontracks/')
source('semitracks.R')

example <- 2

if (example==1){
    fname <- 'track-lengths/fct_c_semi.csv'
    confined <- FALSE
    skip <- 0
    cols <- c(5,9)
} else if (example==2){
    fname <- 'track-lengths/UM10-10 Brogo 3D Semi Tracks with Dper.csv'
    confined <- FALSE
    skip <- 2
    cols <- c(5,9)
} else if (example==3){
    fname <- 'track-lengths/UM10-10 Brogo 3D Confined Lengths.csv'
    confined <- TRUE
    skip <- 2
    cols <- c(5,9)
} else if (example==4){
    fname <- '~/Dropbox/Programming/R/semistijn.csv'
    confined <- FALSE
    skip <- 5
    cols <- c(6,9)
}

# open the data file
dat <- read.data(fname,confined=confined,cutoff=5,skip=skip,cols=cols)
fit <- invert(dat$l,confined=confined,ncomp=3)
mod <- forward(fit$P,fit$M,fit$S)
plotModel(mod,dat)
