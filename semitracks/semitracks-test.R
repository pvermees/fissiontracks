rm(list=ls())
graphics.off()
setwd('~/Documents/Programming/R/fissiontracks/semitracks/')
source('semitracks.R')

example <- 4

if (example==0){
    l <- seq(from=0,to=20,length.out=50)
    M <- 8
    S <- 1.79
    dS <- 1
    fsa <- lapply(l,getfs_l,P=1,M=M,S=S)
    fsb <- lapply(l,getfs_l,P=1,M=M,S=S-dS)
    fla <- lapply(l,getfc_l,P=1,M=M,S=S)
    flb <- lapply(l,getfc_l,P=1,M=M,S=S-dS)
    matplot(cbind(l,l,l,l),cbind(fsa,fla,fsb,flb),
            lty=c(1,1,2,2),col=c('black','blue','black','blue'),
            type='l',xlab='track length (l)',ylab='f(l)')
    legend('topright',legend=c(paste0('semitracks (S=',S,')'),
                               paste0('confined tracks (S=',S,')'),
                               paste0('semitracks (S=',S-dS,')'),
                               paste0('confined tracks (S=',S-dS,')')),
           lty=c(1,1,2,2),col=c('black','blue','black','blue'))
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
