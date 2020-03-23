# fname: the name of a .csv file containing the length measurements
# confined = FALSE if fname contains semi-tracks,
#          = TRUE if it contains confined tracks
# header: does the .csv file contain a header?
# skip: how many lines does the header contain
# cols: column numbers of the true lengths and angles to the c-axis
read.data <- function(fname,confined=FALSE,header=TRUE,skip=0,cols=c(5,9)){
    dat <- read.csv(file=fname,header=header,skip=skip)
    nr <- nrow(dat)
    out <- matrix(0,nr,2)
    colnames(out) <- c('length','angle')
    out[,'length'] <- dat[,cols[1]]
    out[,'angle'] <- dat[,cols[2]]*pi/180
    if (confined) class(out) <- 'confined'
    else class(out) <- 'semi'
    out
}

# dat: input data read by the read.data function
# ncomp: the number of components to fit
invert <- function(dat,confined=FALSE,ncomp=2,mM=5,MM=16){
    # initialise with equal logratios
    pms <- rep(0,2*ncomp)
    if (confined) fit <- optim(pms,cmisfit,dat=dat)
    else fit <- optim(pms,smisfit,dat=dat,mM=mM,MM=MM)
    if (ncomp==1) P <- 1
    else P <- p2P(fit$par[1:(ncomp-1)])
    m <- fit$par[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- fit$par[2*ncomp]
    S <- s2S(s)
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
smisfit <- function(pms,dat,mM=5,MM=16,mS=0.2,MS=2){
    ncomp <- length(pms)/2
    if (ncomp==1) P <- 1
    else P <- p2P(pms[1:(ncomp-1)])
    m <- pms[ncomp:(2*ncomp-1)]
    M <- m2M(m,mM=mM,MM=MM)
    s <- pms[2*ncomp]
    S <- s2S(s,mS=mS,MS=MS)
    out <- 0
    for (i in 1:ncomp){
        f <- getfs_l_phi(dat=dat,M=M[i],S=S)
        out <- out - sum(log(P[i]*f))
    }
    out
}

getfs_l_phi <- function(dat,M,S){
    l <- dat[,'length']
    phi <- dat[,'angle']
    eta <- getEta(M)
    xi <- getXi(M)
    Ms <- eta+4*xi/(3*pi)
    (1-pnorm(q=l,mean=eta+xi*cos(phi),sd=S))*(sin(phi)^2)*4/(pi*Ms)
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
plotModel <- function(m,d=NULL){
    xlab <- expression(paste("length [",mu,"m]"))
    if (is.null(d)){
        par(pty='s',mfrow=c(1,2))
        plot(m$l,m$fc,type='l',xlab=xlab,ylab='')
        title(main='confined tracks')
        plot(m$l,m$fs,type='l',xlab=xlab,ylab='')
        title(main='semi-tracks')
    } else {
        h <- hist(d$l)
        par(pty='s',mfrow=c(2,2))
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
# S = vector with the standard deviations of the length peaks
# P, M and S must have the same lengths
forward <- function(P,M,S,mM=5,MM=16) {
    # calculate the confined and semi-track length distributions
    l <- seq(from=0,to=20,by=0.01) # confined fission track lengths to evaluate
    fc <- l
    fs <- l
    for (i in 1:length(l)){       # loop through all the lengths
        fc[i] <- getfc_l(l[i],P,M,S)  # calculate the pdf of the confined track lengths
        fs[i] <- getfs_l(l[i],P,M,S)  # calculate the pdf of the semi-track lengths
    }
    pms <- c(P2p(P),M2m(M,mM,MM),S2s(S))
    list(l=l,fs=fs,fc=fc,P=P,M=M,S=S,pms=pms)
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
S2s <- function(S,mS=0.2,MS=2){
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
s2S <- function(s,mS=0.2,MS=2){
    mS + (MS-mS)*exp(s)/(1+exp(s))
}
