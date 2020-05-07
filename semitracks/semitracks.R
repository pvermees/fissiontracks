# fname: the name of a .csv file containing the length measurements
# confined = FALSE if fname contains semi-tracks,
#          = TRUE if it contains confined tracks
# header: does the .csv file contain a header?
# skip: how many lines does the header contain
# cols: column numbers of the true lengths and angles to the c-axis
read.data <- function(fname,confined=FALSE,header=TRUE,skip=0,cols=c(5,9)){
    dat <- read.csv(file=fname,header=header,skip=skip)
    l <- dat[,cols[1]]
    a <- dat[,cols[2]]*pi/180
    out <- cbind(l,a)
    colnames(out) <- c('length','angle')
    if (confined) class(out) <- 'confined'
    else class(out) <- 'semi'
    out
}

# dat: input data read by the read.data function
# ncomp: the number of components to fit
# mM and MM: minimum and maximum limit for the modal lengths
# mS and MS: minimum and maximum limit for the peak widths
# r0: cutoff value below which the semi-track lengths are deemed unreliable
invert <- function(dat,confined=FALSE,ncomp=2,mM=5,MM=16,mS=0.2,MS=2.0,r0=0){
    # initialise with equal logratios
    pms <- rep(0,2*ncomp)
    if (confined) fit <- optim(pms,cmisfit,dat=dat)
    else fit <- optim(pms,smisfit,dat=dat,mM=mM,MM=MM,mS=mS,MS=MS,r0=r0)
    if (ncomp==1) P <- 1
    else P <- p2P(fit$par[1:(ncomp-1)])
    m <- fit$par[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- fit$par[2*ncomp]
    S <- s2S(s,mS=mS,MS=MS)
    list(P=P,M=M,S=S,r0=r0)
}

# maximum (negative) likelihood misfit function for confined tracks
cmisfit <- function(pms,dat,mM=5,MM=16,mS=0.2,MS=2.0){
    ncomp <- length(pms)/2
    if (ncomp==1) P <- 1
    else P <- p2P(pms[1:(ncomp-1)])
    m <- pms[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- pms[2*ncomp]
    S <- s2S(s,mS=mS,MS=MS)
    -sum(log(getfc_l_phi(l=dat[,'length'],phi=dat[,'angle'],P=P,M=M,S=S)))
}

# maximum (negative) likelihood misfit function for semi-tracks
# ^f
# |  ---------
# | |     |   \            |
# | |     |    \           |
# | |     |     \          |
# | |     |      ----\     |
# | |     |           \    |
# | |     |            \   |
# --r0---mM----------------MM-----> l
smisfit <- function(pms,dat,mM=5,MM=16,mS=0.9,MS=1.1,r0=0){
    longenough <- dat[,'length']>r0
    l <- dat[longenough,'length']
    a <- dat[longenough,'angle']
    ncomp <- length(pms)/2
    if (ncomp==1) P <- 1
    else P <- p2P(pms[1:(ncomp-1)])
    m <- pms[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- pms[2*ncomp]
    S <- s2S(s,mS=mS,MS=MS)
    Fs_r0 <- getFs_r0(P=P,M=M,S=S,r0=r0)
    -sum(log(getfs_l_phi(l=l,phi=a,P=P,M=M,S=S)/Fs_r0))
}

# P = vector with the proportions of the length peaks
# M = vector with the means of the length peaks
# S = the standard deviation of the length peaks
# P and M must have the same lengths
getfs_l_phi <- function(l,phi,P,M,S){
    out <- 0
    for (i in 1:length(P)){
        eta <- getEta(M[i])
        xi <- getXi(M[i])
        Ms <- eta+4*xi/(3*pi)
        mu <- eta+xi*cos(phi)
        out <- out + P[i]*(1-pnorm(q=l,mean=mu,sd=S))*(sin(phi)^2)*4/(pi*Ms)
    }
    out
}

getFs_r0 <- function(P,M,S,r0=0){
    integrand <- Vectorize(
        function(r0,phi=phi,P=P,M=M,S=S) {
            integrate(getfs_l_phi,lower=0,upper=r0,phi=phi,P=P,M=M,S=S)$value
        }, vectorize.args='phi')
    1 - integrate(integrand,lower=0,upper=pi/2,r0=r0,P=P,M=M,S=S)$value
}

# equation 7.32
# probability density at semi-track length tt
getfs_l <- function(tt,P,M,S){
    integrand <- function(phi,tt,P,M,S) { getfs_l_phi(l=tt,phi=phi,P=P,M=M,S=S) }
    integrate(integrand,lower=0,upper=pi/2,tt=tt,P=P,M=M,S=S)$value
}

# equation 7.51
# get the probability density of length l and angle(s) phi
# M, S, l, theta and f are scalars
getfc_l_phi <- function(l,phi,P,M,S){
    out <- 0
    for (i in 1:length(P)){
        eta <- getEta(M[i])
        xi <- getXi(M[i])
        out <- out + P[i]*dnorm(l,mean=eta+xi*cos(phi),sd=S)
    }
    out
}

# equation 2.4
# get the probability density at confined FT length l integrating over phi
# l is a scalar, while P, M and S are vectors of equal length
getfc_l <- function(l,P,M,S){
    integrand <- function(phi,l,P,M,S) {sin(phi)*getfc_l_phi(l,phi,P,M,S)}
    integrate(integrand,lower=0,upper=pi/2,l=l,P=P,M=M,S=S)$value
}

# equation 7.53
getEta <- function(mu){
    -5.18 + 1.42*mu - sqrt(11.48 - 2.13*mu + 0.1*mu^2)
}

# equation 7.53
getXi <- function(mu){
    2*(mu-getEta(mu))
}

# middle of page 140
getMs <- function(mu){
    getEta(mu)+4*getXi(mu)/(3*pi)
}

getELA <- function(mu,phi){
    getEta(mu) + getXi(mu)*cos(phi)
}

# fit = forward model, dat = data
plotModel <- function(fit,P=1,M=10,S=1,r0=0,dat=NULL){
    if (missing(fit)){
        fit <- list(P=P,M=M,S=S,r0=r0)
    }
    par(pty='s',mfrow=c(2,2))
    m <- forward(fit)
    xlab <- expression(paste("length [",mu,"m]"))
    hasdat <- !is.null(dat)
    confined <- hasdat && is(dat,'confined')
    semi <- hasdat && is(dat,'semi')
    if (hasdat){
        longenough <- dat[,'length']>fit$r0
        d <- dat[longenough,]
        h <- hist(d[,'length'],plot=FALSE)
    }
    if (confined){
        plot(h,freq=FALSE,main='',ylim=range(h$density,m$fc_l),xlab=xlab)
        lines(m$l,m$fc_l,type='l',xlab=xlab,ylab='')
    } else {
        plot(m$l,m$fc_l,type='l',xlab=xlab,ylab='')
    }
    title(main='confined tracks')
    if (semi){
        plot(h,freq=FALSE,main='',ylim=range(h$density,m$fs_l),xlab=xlab)
        lines(m$l,m$fs_l,type='l',xlab=xlab,ylab='')
    } else {
        plot(m$l,m$fs_l,type='l',xlab=xlab,ylab='')
    }
    title(main='semi-tracks')
    if (confined){
        plot(d[,'length'],d[,'angle']*180/pi,pch=21,col='grey',
             xlab=xlab,ylab='angle (degrees)')
        contour(x=m$l,y=m$phi*180/pi,z=m$fc_l_phi,
                xlab=xlab,ylab='angle (degrees)',add=TRUE)

    } else {
        contour(x=m$l,y=m$phi*180/pi,z=m$fc_l_phi,
                xlab=xlab,ylab='angle (degrees)',add=FALSE)
    }
    title(main='confined tracks')
    if (semi){
        plot(d[,'length'],d[,'angle']*180/pi,pch=21,col='grey',
             xlab=xlab,ylab='angle (degrees)')
        contour(x=m$l,y=m$phi*180/pi,z=m$fs_l_phi,add=TRUE)
    } else {
        contour(x=m$l,y=m$phi*180/pi,z=m$fs_l_phi,
                xlab=xlab,ylab='angle (degrees)')
    }
    title(main='semi-tracks')
}

# nn = number of points at which to evaluate fc and fs
forward <- function(fit,nn=50) {
    P <- fit$P
    M <- fit$M
    S <- fit$S
    r0 <- fit$r0
    nn <- 100
    out <- list()
    out$l <- seq(from=0,to=20,length.out=nn)
    out$phi <- seq(from=0,to=pi/2,length.out=nn)
    out$fc_l <- rep(0,nn)
    out$fs_l <- rep(0,nn)
    out$fc_l_phi <- matrix(0,nn,nn)
    out$fs_l_phi <- matrix(0,nn,nn)
    out$Fs_r0 <- getFs_r0(P=P,M=M,S=S,r0=r0)
    for (i in 1:nn){
        for (j in 1:nn){
            out$fc_l_phi[i,j] <- getfc_l_phi(l=out$l[i],phi=out$phi[j],P=P,M=M,S=S)
            out$fs_l_phi[i,j] <- getfs_l_phi(l=out$l[i],phi=out$phi[j],P=P,M=M,S=S)/out$Fs_r0
        }
        out$fc_l[i] <- getfc_l(l=out$l[i],P=P,M=M,S=S)
        out$fs_l[i] <- getfs_l(tt=out$l[i],P=P,M=M,S=S)/out$Fs_r0
    }
    out
}

## log(ratio) transformations:
P2p <- function(P){
    nn <- length(P)
    log(P[1:(nn-1)]) - log(P[nn])
}
M2m <- function(M,mM=5,MM=16){
    nM <- length(M)
    dM <- diff(c(mM,M))
    log(dM) - log(MM-M[nM])
}
S2s <- function(S,mS=0.9,MS=1.1){
    log(S-mS) - log(MS-S)
}

## inverse log(ratio) transformations:
p2P <- function(p){
    c(exp(p),1)/(sum(exp(p))+1)
}
# m = log of the ratios of the differences between the peak locations
# mM, MM = minimum and maximum track length
m2M <- function(m,mM=5,MM=16){
    dM <- (MM-mM)*exp(m)/(1+sum(exp(m)))
    mM + cumsum(dM)
}
s2S <- function(s,mS=0.9,MS=1.1){
    mS + (MS-mS)*exp(s)/(1+exp(s))
}

# create a synthetic dataset of horizontally confined fission track lengths
hcft <- function(fit,nn=100){
    misfit <- function(l,q,fit){
        Q <- 0
        for (i in 1:length(fit$P)){
            Q <- Q + fit$P[i]*pnorm(q=l,mean=fit$M[i],sd=fit$S)
        }
        (Q - q)^2
    }
    out <- rep(0,nn)
    quantile <- seq(from=0,to=1,length.out=(nn+2))[2:(nn+1)]
    for (i in 1:nn){
        out[i] <- optimise(misfit,interval=c(0,20),q=quantile[i],fit=fit)$minimum
    }
    out
}
