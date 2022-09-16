rm(list=ls())
graphics.off()
setwd('~/Documents/Programming/R/fissiontracks/semitracks/')
source('semitracks.R')

example <- 0

if (example==0){
    i <- 2
    M <- c(15,8,11.5)
    S <- c(0.87,1.79,0.97)
    l <- seq(from=0,to=20,length.out=200)
    fs <- lapply(l,getfs_l,P=1,M=M[i],S=S[i])
    fl <- lapply(l,getfc_l,P=1,M=M[i],S=S[i])
    matplot(cbind(l,l),cbind(fs,fl),lty=1,col='black',type='l')
} else if (example==1){
    fname <- 'track-lengths/fct_c_semi.csv'
    confined <- FALSE
    skip <- 0
    cols <- c(5,8)
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
    fname <- 'track-lengths/GV15JAW261_confined-track_lengths.csv'
    confined <- TRUE
    skip <- 0
    cols <- c(5,9)

    fit <- list(P=c(0.7424027,0.2575973),
                M=c(10.65074,12.98421),
                S=0.9704537,
                r0=1)

} else if (example==5){
    fname <- 'track-lengths/GV15JAW261_semi-track_lengths.csv'
    confined <- FALSE
    skip <- 0
    cols <- c(5,9)

    fit <- list(P=c(0.2088414,0.7911586),
                M=c(11.56809,13.25569),
                S=0.7923793,
                r0=1)
    
}


if (example>0){
    dat <- read.data(fname,confined=confined,skip=skip,cols=cols)
    if (!exists('fit')){
        fit <- invert(dat,confined=confined,ncomp=2,r0=1)
    }
    plotModel(fit=fit,dat=dat)
    dat <- hcft(fit=fit,nn=200)
}
