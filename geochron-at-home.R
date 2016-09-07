rm(list=ls())

library('jpeg')

# ij = c(grain number, image number)
load.grain <- function(dl,ij=c(1,1),mica=FALSE,vertices=NULL,new=FALSE){
    if (!new){
        par(mar=c(1,1,1,1),xpd=NA,new=new)
        plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab='',ylab='')
        text(0.75,0,'new',pos=1)
        text(0.25,0,'delete',pos=1)
        text(0.75,1,'mica',pos=3)
        text(0.25,1,'grain',pos=3)
        text(0,0.75,'next',pos=2)
        text(0,0.25,'prev',pos=2)
        text(1,0.75,'up',pos=4)
        text(1,0.25,'down',pos=4)
        text(1,-0.02,'exit',pos=4)
        title(main=paste0('Grain ',ij[1]))
    }
    i <- ij[1] # grain number
    j <- ij[2] # image number
    roi <- vertices[[i]]
    ng <- length(dl) # number of grains
    if (i<1) i <- 1
    if (i>ng) i <- ng
    if (mica) fl <- list.files(dl[i],pattern="^MicaStack-")
    else fl <- list.files(dl[i],pattern="^Stack-")
    fl <- fl[sort.fl(fl,mica)]
    if (mica) fl <- c("MicaReflStackFlat.jpeg",fl)
    else fl <- c("ReflStackFlat.jpeg",fl)
    ni <- length(fl) # number of images
    if (j<1) j <- 1
    if (j>ni) j <- ni
    f <- paste(dl[i],fl[j],sep='/') # file path
    img <- readJPEG(f)
    height <- dim(img)[1]
    width <- dim(img)[2]
    if (mica) rasterImage(img[,width:1,],0,0,1,1)
    else rasterImage(img,0,0,1,1)
    draw.roi(roi)
    c(i,j,height,width)
}

sort.fl <- function(fl,mica=FALSE){
    nf <- length(fl)
    ii <- rep(0,nf)
    if (mica) start <- 10
    else start <- 7
    for (i in 1:nf){
        nn <- nchar(fl[i])
        ii[i] <- as.numeric(substr(fl[i],start=start,stop=nn-5))
    }
    sort(ii,index.return=TRUE,decreasing=mica)$ix
}

init.roi <- function(dl){
    ng <- length(dl)
    out <- list()
    for (i in 1:ng){
        fl <- list.files(dl[i])
        roi <- list(x=NULL,y=NULL)
        if ("rois.json" %in% fl) {
            fname <- paste0(dl[i],"/rois.json")
            json <- paste(readLines( fname, warn=FALSE ),collapse="")
            w <- as.numeric(unlist(strsplit(json,'"image_width":|,'))[2])
            h <- as.numeric(unlist(strsplit(json,'"image_height":|,'))[3])
            vertstr <- unlist(strsplit(json,'\\[\\[|\\]\\]'))[2]
            corners <- unlist(strsplit(vertstr,'\\],\\['))
            for (j in 1:length(corners)){
                coords <- as.numeric(unlist(strsplit(corners[j],',')))
                if (!any(is.na(coords))){
                    roi$x <- c(roi$x,coords[1]/w)
                    roi$y <- c(roi$y,coords[2]/h)
                }
            }
        }
        out[[i]] <- roi
    }
    out
}

donext <- function(xy){
    (xy$x<0) & (xy$y>=0.5) & (xy$y<1)
}
doprev <- function(xy){
    (xy$x<0) & (xy$y<0.5) & (xy$y>0)
}
dograin <- function(xy){
    (xy$x<=0.5) & (xy$x>0) & (xy$y>1)
}
domica <- function(xy){
    (xy$x>0.5) & (xy$x<1) & (xy$y>1)
}
doup <- function(xy){
    (xy$x>1) & (xy$y>0.5) & (xy$y<1)
}
dodown <- function(xy){
    (xy$x>1) & (xy$y<=0.5) & (xy$y>0)
}
doexit <- function(xy){
    (xy$x>1) & (xy$y<0)
}
donew <- function(xy){
    (xy$x>0.5) & (xy$x<1) & (xy$y<0)
}
dodelete <- function(xy){
    (xy$x<=0.5) & (xy$x>0) & (xy$y<0)
}

geochron.at.home <- function(dir){
    mica <- FALSE
    dl <- list.dirs(dir,recursive=FALSE) # directory list
    vertices <- init.roi(dl)
    ij <- load.grain(dl,mica=mica,vertices=vertices,new=FALSE)
    while(TRUE){
        message("Choose one of the following options:\n",
                "b - browse the images\n",
                "r - mark regions of interest\n",
                "s - determine the shift between grain and mica\n",
                "c - count fission tracks\n",
                "x - exit the program and save the output")
        response <- readline()
        if (identical(response,'r') |
            identical(response,'b') |
            identical(response,'c')){
            while(TRUE){
                xy <- locator(1)
                if (dodown(xy)){
                    ij[2] <- ij[2]+1
                    ij <- load.grain(dl,ij,mica,vertices,new=TRUE)
                } else if (doup(xy)){
                    ij[2] <- ij[2]-1
                    ij <- load.grain(dl,ij,mica,vertices,new=TRUE)
                } else if (doprev(xy) & ij[1]>1){
                    if (!identical(response,'b'))
                        save.to.json(dl[ij[1]],ij[3],ij[4],vertices[[ij[1]]])
                    ij[1:2] <- c(ij[1]-1,1)
                    ij <- load.grain(dl,ij,mica,vertices,new=FALSE)
                } else if (donext(xy) & ij[1]<length(vertices)){
                    if (!identical(response,'b'))
                        save.to.json(dl[ij[1]],ij[3],ij[4],vertices[[ij[1]]])
                    ij[1:2] <- c(ij[1]+1,1)
                    ij <- load.grain(dl,ij,mica,vertices,new=FALSE)
                } else if (donew(xy)){
                    # TODO
                } else if (doexit(xy)){
                    if (!identical(response,'b'))
                        save.to.json(dl[ij[1]],ij[3],ij[4],vertices[[ij[1]]])
                    break;
                } else if (dodelete(xy)){
                    roi <- vertices[[ij[1]]]
                    nc <- length(roi$x)
                    if (nc==1){
                        roi$x <- NULL
                        roi$y <- NULL
                    } else {
                        roi$x <- roi$x[1:(nc-1)]
                        roi$y <- roi$y[1:(nc-1)]
                    }
                    vertices[[ij[1]]] <- roi
                    load.grain(dl,ij,mica,vertices,new=TRUE)
                } else if (domica(xy)){
                    mica <- TRUE
                    ij[2] <- 1
                    load.grain(dl,ij,mica,vertices,new=FALSE)
                } else if (dograin(xy)){
                    mica <- FALSE
                    ij[2] <- 1
                    load.grain(dl,ij,mica,vertices,new=FALSE)            
                } else if (!identical(response,'b')) {
                    vertices[[ij[1]]]$x <- c(vertices[[ij[1]]]$x,xy$x)
                    vertices[[ij[1]]]$y <- c(vertices[[ij[1]]]$y,1-xy$y)
                    load.grain(dl,ij,mica,vertices,new=TRUE)
                }
            }
        } else if (identical(response,'x')){
            graphics.off()
            break
        }
    }
}

draw.roi <- function(roi){
    nc <- length(roi$x)
    if (nc>0){
        for (i in 2:nc){
            lines(roi$x[(i-1):i],1-roi$y[(i-1):i],col='white')
        }
        lines(roi$x[c(1,nc)],1-roi$y[c(1,nc)],col='white')
    }
}

save.to.json <- function(dir,width,height,roi){
    out <- "{"
    out <- paste0(out,'"image_width":',width,',')
    out <- paste0(out,'"image_height":',height,',')
    out <- paste0(out,'"regions":[')
    nc <- length(roi$x) # number of corners
    out <- paste0(out,'{"vertices":[')
    for (j in 1:nc){
        out <- paste0(out,'[',round(width*roi$x[j]),
                      ',',round(height*roi$y[j]),']')
        if (j<nc) out <- paste0(out,",")
    }
    out <- paste0(out,'],"shift":[0,0]}]}')
    sink(paste0(dir,'/rois.json'))
    cat(out)
    sink()
}

geochron.at.home("LU324-2-DUR")
geochron.at.home("LU324-6-FCT")
