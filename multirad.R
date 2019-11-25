library(IsoplotR)

FTpool <- function(...){
    list(...)
}

# bp = beta-phi, x = pooled FT data
# i = peak number, j = irradiation number, u = grain number
fiju <- function(bi,x,j,u){
    stats::dbinom(yju(x,j,u),mju(x,j,u),thetaij(bi,x,j))
}

aiu <- function(b,x,i,j,u){
    yju(x,j,u) - thetaij(b[i],x,j)*mju(x,j,u)
}

biu <- function(b,x,i,j,u){
    aiu(b,x,i,j,u)^2 - thetaij(b[i],x,j) *
                       (1-thetaij(b[i],x,j))*mju(x,j,u)
}

get.piju <- function(b,p,x,i,j,u){
    den <- 0
    for (k in 1:nd(x)){
        den <- den + p[k]*fiju(b[k],x,j,u)
    }
    num <- p[i]*fiju(b[i],x,j,u)
    num/den
}

get.piu <- function(b,p,x){
    k <- length(b)
    piu <- NULL
    for (j in 1:nd(x)){
        for (u in 1:nj(x,j)){
            newrow <- rep(0,k)
            for (i in 1:k){
                newrow[i] <- get.piju(b,p,x,i,j,u)
            }
            piu <- rbind(piu,newrow)
        }
    }
    piu
}

mju <- function(x,j,u){
    sum(x[[j]]$x[u,])
}

yju <- function(x,j,u){
    x[[j]]$x[u,'Ns']
}

update.p <- function(x,b,p,piu){
    out <- 0*p
    k <- length(b)
    for (i in 1:k){
        rn <- 0 # row number in piu matrix
        for (j in 1:nd(x)){
            for (u in 1:nj(x,j)){
                rn <- rn + 1
                out[i] <- out[i] + piu[rn,i]
            }
        }
    }
    out/rn
}

pii <- function(x,b,p,i){
    out <- 0
    n <- 0
    for (j in 1:nd(x)){
        for (u in 1:nj(x,j)){
            n <- n + nj[j]
            out <- out + get.piju(b,p,x,i,j,u)
        }
    }
    out
}

beta.misfit <- function(bi,x,b,p,piu,i){
    LHS <- 0
    RHS <- 0
    rn <- 0
    for (j in 1:nd(x)){
        for (u in 1:nj(x,j)){
            rn <- rn + 1
            LHS <- LHS + piu[rn,i]*yju(x,j,u)
            RHS <- RHS + thetaij(bi,x,j)*piu[rn,i]*mju(x,j,u)
        }
    }
    ((LHS-RHS)/(LHS+RHS))^2
}

update.b <- function(x,b,p,piu){
    k <- length(b)
    out <- b
    for (i in 1:k){
        out[i] <- optimise(f=beta.misfit,interval=c(-30,30),
                           x=x,b=b,p=p,piu=piu,i=i)$minimum
    }
    out
}

nd <- function(x){
    length(x)
}

nj <- function(x,j){
    nrow(x[[j]]$x)
}

nu <- function(x){
    out <- 0
    for (j in 1:nd(x)){
        out <- out + nj(x,j)
    }
    out
}

thetaij <- function(bi,x,j){
    exp(bi)/(exp(bi)+x[[j]]$rhoD[1])
}

peakfit <- function(x,k=2){
    b <- init.b(x,k=k)
    p <- rep(1,k)/k
    for (iter in 1:100){
        piu <- get.piu(b,p,x)
        p <- update.p(x,b,p,piu)
        b <- update.b(x,b,p,piu)
    }
    l238 <- 0.000155125
    zeta <- x[[1]]$zeta[1]/1e6
    age <- log(1 + 0.5*l238*zeta*exp(b))/l238
    out <- rbind(age,p)
    rownames(out) <- c('t','p')
    out
}

init.b <- function(x,k=2){
    b <- rep(0,nu(x))
    rn <- 0
    for (i in 1:nd(x)){
        rhoD <- x[[i]]$rhoD[1]
        for (u in 1:nj(x,i)){
            rn <- rn + 1
            NsNi <- (x[[i]]$x[u,'Ns']+0.5)/(x[[i]]$x[u,'Ni']+0.5)
            b[rn] <- log(rhoD*NsNi)
        }
    }
    seq(min(b),max(b),length.out=k)
}
