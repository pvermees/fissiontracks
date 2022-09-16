# equation 7.32
# returns the probability density at semi-track length tt
# P = vector with the proportions of the length peaks
# M = vector with the means of the length peaks
# S = the standard deviation of the length peaks
# P and M must have the same lengths
getfs_l <- function(tt,P,M,S){
    integrand <- function(phi,tt,P,M,S) { getfs_l_phi(l=tt,phi=phi,P=P,M=M,S=S) }
    integrate(integrand,lower=0,upper=pi/2,tt=tt,P=P,M=M,S=S)$value
}

# helper function for getfs_l
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

# equation 2.4
# get the probability density at confined FT length l integrating over phi
# l is a scalar, while P, M and S are vectors of equal length
getfc_l <- function(l,P,M,S){
    integrand <- function(phi,l,P,M,S) {sin(phi)*getfc_l_phi(l,phi,P,M,S)}
    integrate(integrand,lower=0,upper=pi/2,l=l,P=P,M=M,S=S)$value
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

# equation 7.53
getEta <- function(mu){
    -5.18 + 1.42*mu - sqrt(11.48 - 2.13*mu + 0.1*mu^2)
}

# equation 7.53
getXi <- function(mu){
    2*(mu-getEta(mu))
}

# an attempt to reproduce the first two panels of Figure 7.5.b:
fs <- lapply(l,getfs_l,P=1,M=8,S=1.79)
fl <- lapply(l,getfc_l,P=1,M=8,S=1.79)
matplot(cbind(l,l),cbind(fs,fl),
        lty=1,col='black',type='l',xlab='track length (l)',ylab='f(l)')
