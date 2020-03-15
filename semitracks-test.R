rm(list=ls())
graphics.off()
setwd('~/Documents/Programming/R/fissiontracks/')
source('semitracks.R')

example <- 4

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
    fname <- '~/Desktop/Length_data_semi.csv'
    confined <- FALSE
    skip <- 5
    cols <- c(6,9)
}

dat <- read.data(fname,confined=confined,skip=skip,cols=cols)

option <- 2

if (option==1) {
    P <- c(0.3,0.7)
    M <- c(8,13)
    S <- 1
    fm <- forward(P,M,S)
    plotModel(fm)
    l <- simulate(fm,n=1000,confined=FALSE)
    fit <- invert(l,confined=FALSE,ncomp=2)
} else if (option==2){
    fit <- invert(dat$l,confined=confined,ncomp=1)
    fm <- forward(fit$P,fit$M,fit$S)
    plotModel(fm,dat)
} else if (option==3){ # plot fs against l for different phi
    nl <- 50
    np <- 10
    l <- seq(from=0,to=15,length.out=nl)
    phi <- seq(from=pi/2,to=0,length.out=np)
    M <- 10
    S <- 1
    for (i in 1:np){
        dat <- cbind(l,phi[i])
        colnames(dat) <- c('length','angle')
        f <- getfsLA(dat,M,S)
        if (i==1) plot(l,f,type='l')
        else lines(l,f,type='l')
    }
} else if (option==4){ # plot E[l,mu] on angle plot
    nm <- 3
    M <- seq(from=15,to=7,length.out=nm)
    np <- 100
    phi <- seq(from=0,to=pi/2,length.out=np)
    plot(c(0,20),c(0,100),type='n',xlab='length',ylab='angle')
    for (i in 1:nm){
        l <- getElphi(mu=M[i],phi)
        lines(l,phi*180/pi,type='l')
    }
} else if (option==5){ # plot fs against l for different M
    nl <- 50
    nm <- 3
    l <- seq(from=0,to=20,length.out=nl)
    phi <- pi/2
    M <- seq(from=7,to=15,length.out=nm)
    S <- 1
    for (i in 1:nm){
        dat <- cbind(l,phi)
        colnames(dat) <- c('length','angle')
        f <- getfsLA(dat,M[i],S)
        if (i==1) plot(l,f,type='l')
        else lines(l,f,type='l')
    }
} else if (option==6){
    fit <- optim(log(c(10,1)),smisfit,dat=dat)
}

#mod <- forward(fit$P,fit$M,fit$S)
#plotModel(mod,dat)
