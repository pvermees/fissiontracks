# fname: the name of a .csv file containing the length measurements
# confined = FALSE if fname contains semi-tracks,
#          = TRUE if it contains confined tracks
# header: does the .csv file contain a header?
# skip: how many lines does the header contain
# cols: column numbers of the true lengths and angles to the c-axis
read.data <- function(fname,confined=FALSE,header=TRUE,skip=0,cols=c(5,9),cutoff=0){
    dat <- read.csv(file=fname,header=header,skip=skip)
    l <- dat[,cols[1]]
    a <- dat[,cols[2]]*pi/180
    long <- (l>cutoff)
    out <- cbind(l[long],a[long])
    colnames(out) <- c('length','angle')
    if (confined) class(out) <- 'confined'
    else class(out) <- 'semi'
    out
}

# dat: input data read by the read.data function
# ncomp: the number of components to fit
invert <- function(dat,confined=FALSE,ncomp=2,mM=5,MM=16,mS=0.2,MS=2.0){
    # initialise with equal logratios
    pms <- rep(0,2*ncomp)
    if (confined) fit <- optim(pms,cmisfit,dat=dat)
    else fit <- optim(pms,smisfit,dat=dat,mM=mM,MM=MM,mS=mS,MS=MS)
    if (ncomp==1) P <- 1
    else P <- p2P(fit$par[1:(ncomp-1)])
    m <- fit$par[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- fit$par[2*ncomp]
    S <- s2S(s,mS=mS,MS=MS)
    list(P=P,M=M,S=S)
}

# maximum (negative) likelihood misfit function for semi-tracks
# ^f
# |----------
# |      |   \            |
# |      |    \           |
# |      |     \          |
# |      |      ----\     |
# |      |           \    |
# |      |            \   |
# ------mM----------------MM-----> l
smisfit <- function(pms,dat,mM=5,MM=16,mS=0.9,MS=1.1){
    ncomp <- length(pms)/2
    if (ncomp==1) P <- 1
    else P <- p2P(pms[1:(ncomp-1)])
    m <- pms[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- pms[2*ncomp]
    S <- s2S(s,mS=mS,MS=MS)
    out <- 0
    for (i in 1:ncomp){
        f <- getfs_l_phi(l=dat[i,'length'],phi=dat[i,'angle'],M=M[i],S=S)
        out <- out - log(sum(P[i]*f))
    }
    out
}

getfs_l_phi <- function(l,phi,P,M,S){
    out <- 0
    for (i in 1:length(P)){
        eta <- getEta(M[i])
        xi <- getXi(M[i])
        Ms <- eta+4*xi/(3*pi)
        out <- out + P[i]*(1-pnorm(q=l,mean=eta+xi*cos(phi),sd=S))*(sin(phi)^2)*4/(pi*Ms)
    }
    out
}

# equation 7.32
# probability density at semi-track length tt
getfs_l <- function(tt,P,M,S){
    out <- 0 # initialise
    for (i in 1:length(P)){
        integrand <- function(phi,tt,Pi,Mi,S) { getfs_l_phi(l=tt,phi=phi,Pi,Mi,S) }
        out <- out + integrate(integrand,lower=0,upper=pi/2,
                               tt=tt,Pi=P[i],Mi=M[i],S=S)$value
    }
    out
}

# equation 7.51
# get the probability density of length l and angle(s) phi
# M, S, l, theta and f are scalars
getfc_l_phi <- function(l,phi,P,M,S){
    out <- 0
    for (i in length(P)){
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
    out <- 0
    for (i in 1:length(P)){
        integrand <- function(phi,l,Pi,Mi,S) {sin(phi)*getfc_l_phi(l,phi,P[i],M[i],S)}
        out <- out + integrate(integrand,lower=0,upper=pi/2,
                               l=l,Pi=P[i],Mi=M[i],S=S)$value
    }
    out
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

# m = forward model, d = data
plotModel <- function(P,M,S,d=NULL){
    par(pty='s',mfrow=c(2,2))
    m <- forward(P=P,M=M,S=S)
    xlab <- expression(paste("length [",mu,"m]"))
    if (is.null(d)){
        plot(m$l,m$fc_l,type='l',xlab=xlab,ylab='')
        title(main='confined tracks')
        plot(m$l,m$fs_l,type='l',xlab=xlab,ylab='')
        title(main='semi-tracks')
    } else {
        h <- hist(d$l)
        plot(d$l,d$a,type='p', xlab=xlab, ylab='angle to C')
        plot(h,freq=FALSE,main='',ylim=range(h$density,m$fs),xlab=xlab)
        if (d$confined){
            lines(m$l,m$fc)
            title(smisfit(m$pms,d$l))
            plot(m$l,m$fs,type='l',xlab=xlab)
        } else {
            cutoff <- min(d$l)
            fsL.of.uncounted <- getfs_l(cutoff,P=m$P,M=m$M,S=m$S)
            corr.factor <- 1/(1-cutoff*fsL.of.uncounted)
            lines(m$l,m$fs*corr.factor)
            title(smisfit(m$pms,d$l))
            plot(m$l,m$fc,type='l',xlab=xlab,ylab='Density')
        }
    }
}

# P = vector with the proportions of the length peaks
# M = vector with the means of the length peaks
# S = the standard deviation of the length peaks
# P and M must have the same lengths
# nn = number of points at which to evaluate fc and fs
forward <- function(P,M,S,nn=50) {
    nn <- 100
    out <- list()
    out$l <- seq(from=0,to=20,length.out=nn)
    out$phi <- seq(from=0,to=pi/2,length.out=nn)
    out$fc_l <- rep(0,nn)
    out$fs_l <- rep(0,nn)
    out$fc_l_phi <- matrix(0,nn,nn)
    out$fs_l_phi <- matrix(0,nn,nn)
    for (i in 1:nn){
        out$fc_l[i] <- getfc_l(out$l[i],P=P,M=M,S=S)
        out$fs_l[i] <- getfs_l(out$l[i],P=P,M=M,S=S)
        for (j in 1:nn){
            out$fc_l_phi[i,j] <- getfc_l_phi(l=out$l[i],phi=out$phi[j],P=P,M=M,S=S)
            out$fs_l_phi[i,j] <- getfs_l_phi(l=out$l[i],phi=out$phi[j],P=P,M=M,S=S)
        }
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
