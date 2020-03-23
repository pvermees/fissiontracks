rm(list=ls())
graphics.off()
setwd('~/Documents/Programming/R/fissiontracks/')
source('semitracks2.R')

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
    fname <- '~/Desktop/Length_data_semi.csv'
    confined <- FALSE
    skip <- 5
    cols <- c(6,9)
}

dat <- read.data(fname,confined=confined,skip=skip,cols=cols)

option <- 10

if (option==10){
    nc <- 5
    nr <- 50
    M <- seq(from=5,to=16,length.out=nc)
    S <- seq(from=0.2,to=2,length.out=nc)
    l <- seq(from=0,to=16,length.out=nr)
    phi <- rep(45*pi/180,nr)
    dat <- cbind(l,phi)
    colnames(dat) <- c('length','angle')
    lf <- matrix(0,nr,nc)
    x <- lf
    for (i in 1:nc){
        x[,i] <- l
        lf[,i] <- log(getfs_l_phi(dat,M=10,S=S[i]))
    }
    matplot(x=l,y=lf,type='l')
} else if (option==9){
    fs <- log(getfs_l_phi(dat,M=6,S=10))
    mfs <- min(fs)
    Mfs <- max(fs)
    r <- (fs-mfs)/(Mfs-mfs)
    g <- (Mfs-fs)/(Mfs-mfs)
    b <- 0
    plot(dat[,'length'],dat[,'angle'],pch=19,col=rgb(cbind(r,g,b)))
    title(main=c(mfs,Mfs,sum(fs)))
} else if (option==8){
    nn <- 50
    M <- seq(from=5,to=16,length.out=nn)
    S <- seq(from=0.2,to=2,length.out=nn)
    mf <- matrix(0,nn,nn)
    for (i in 1:nn){
        m <- M2m(M[i])
        for (j in 1:nn){
            s <- S2s(S[i])
            mf[i,j] <- sum(log(getfs_l_phi(dat,M=M[i],S=S[j])))
            #mf[i,j] <- smisfit(c(m,s),dat)
        }
    }
    contour(x=M,y=S,z=mf)
} else if (option==7){
    nn <- 100
    l <- seq(from=0,to=16,length.out=nn)
    a <- seq(from=0,to=pi/2,length.out=nn)
    lf <- matrix(0,nn,nn)
    gdat <- matrix(0,1,2)
    colnames(gdat) <- c('length','angle')
    P1 <- 1
    P2 <- 0.7
    M1 <- 10
    M2 <- 5
    S <- 2
    for (i in 1:nn){
        gdat[1,'length'] <- l[i]
        for (j in 1:nn){
            gdat[1,'angle'] <- a[j]
            lf[i,j] <- log(P1*getfs_l_phi(gdat,M=M1,S=S)) # +
            # log(P2 * getfs_l_phi(gdat,M=M2,S=S))
        }
    }
    pdat <- dat
    pdat[,'angle']<- dat[,'angle']*180/pi
    plot(pdat,pch=21,col='grey')
    contour(x=l,y=a*180/pi,z=lf,
            xlab='length',ylab='angle',add=TRUE)
    p <- 0#P2p(c(P1,P2))
    m <- M2m(M1)#,M2m(c(M1,M2),mM=5,MM=16)
    s <- S2s(S)
    mf <- smisfit(c(m,s),dat,mM=5,MM=16)
    title(main=mf)
} else if (option==6){
    f <- log(getfs_l_phi(dat=dat,M=10,S=1))
    plot(dat[,'angle'],f,type='p')
    title(main=sum(f))
    fit <- invert(dat,confined=confined,ncomp=1,mS=1,MS=1)
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
} else if (option==2){
    fit <- invert(dat$l,confined=confined,ncomp=1)
    fm <- forward(fit$P,fit$M,fit$S)
    plotModel(fm,dat)
} else if (option==1) {
    P <- c(0.3,0.7)
    M <- c(8,13)
    S <- 1
    fm <- forward(P,M,S)
    plotModel(fm)
    l <- simulate(fm,n=1000,confined=FALSE)
    fit <- invert(l,confined=FALSE,ncomp=2)
}

#mod <- forward(fit$P,fit$M,fit$S)
#plotModel(mod,dat)
