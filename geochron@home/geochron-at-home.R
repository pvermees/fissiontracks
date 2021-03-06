rm(list=ls())

library('jpeg')
library('sp')

# ij = c(grain number, image number)
load.grain <- function(dl,ij=c(1,1),mica=FALSE,json=NULL,new=FALSE,
                       lab='Melbourne'){
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
    roi <- json$rois[[i]]
    #print(Polygon(cbind(c(roi$x,roi$x[1]),c(roi$y,roi$y[1])),hole=F)@area)
    spots <- json$tracks[[i]]
    ng <- length(dl) # number of grains
    if (i<1) i <- 1
    if (i>ng) i <- ng
    if (identical(lab,'Melbourne')){
        if (mica) fl <- list.files(dl[i],pattern="^MicaStack-")
        else fl <- list.files(dl[i],pattern="^Stack-")
        fl <- fl[sort.fl(fl,mica)]
        if (mica) fl <- c("MicaReflStackFlat.jpeg",fl)
        else fl <- c("ReflStackFlat.jpeg",fl)
    } else if (identical(lab,'Gent')){
        if (mica){
            fl <- list.files(dl[i],pattern="m12")
        } else {
            fl <- list.files(dl[i],pattern="dia")
            fl <- c(list.files(dl[i],pattern="epi"),fl)
        }
    }
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
    draw.spots(spots)
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

init.json <- function(dl){
    ng <- length(dl)
    out <- list(rois=list(),tracks=list())
    for (i in 1:ng){
        fl <- list.files(dl[i])
        roi <- list(x=NULL,y=NULL)
        tracks <- list(x=NULL,y=NULL)
        if ("rois.json" %in% fl) {
            fname <- paste0(dl[i],"/rois.json")
            json <- paste(readLines( fname, warn=FALSE ),collapse="")
            w <- as.numeric(unlist(strsplit(json,'"image_width":|,'))[2])
            h <- as.numeric(unlist(strsplit(json,'"image_height":|,'))[3])
            vertstr <- unlist(strsplit(json,'\\[\\[|\\]\\]'))[2]
            corners <- unlist(strsplit(vertstr,'\\],\\['))
            trackstr <- unlist(strsplit(json,'\\[\\[|\\]\\]'))[4]
            spots <- unlist(strsplit(trackstr,'\\],\\['))
            for (j in 1:length(corners)){
                coords <- as.numeric(unlist(strsplit(corners[j],',')))
                if (!any(is.na(coords))){
                    roi$x <- c(roi$x,coords[1]/w)
                    roi$y <- c(roi$y,coords[2]/h)
                }
            }
            for (j in 1:length(spots)){
                coords <- as.numeric(unlist(strsplit(spots[j],',')))
                if (!any(is.na(coords))){
                    tracks$x <- c(tracks$x,coords[1]/w)
                    tracks$y <- c(tracks$y,coords[2]/h)
                }
            }
        }
        out$rois[[i]] <- roi
        out$tracks[[i]] <- tracks
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

geochron.at.home <- function(dir,lab='Melbourne'){
    mica <- FALSE
    dl <- list.dirs(dir,recursive=FALSE) # directory list
    json <- init.json(dl)
    ij <- load.grain(dl,mica=mica,json=json,new=FALSE,lab=lab)
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
                    ij <- load.grain(dl,ij,mica,json,new=TRUE,lab=lab)
                } else if (doup(xy)){
                    ij[2] <- ij[2]-1
                    ij <- load.grain(dl,ij,mica,json,new=TRUE,lab=lab)
                } else if (doprev(xy) & ij[1]>1){
                    if (!identical(response,'b'))
                        save.to.json(dl[ij[1]],ij[3],ij[4],
                                     json$rois[[ij[1]]],
                                     json$tracks[[ij[1]]])
                    ij[1:2] <- c(ij[1]-1,1)
                    ij <- load.grain(dl,ij,mica,json,new=FALSE,lab=lab)
                } else if (donext(xy) & ij[1]<length(json$rois)){
                    if (!identical(response,'b'))
                        save.to.json(dl[ij[1]],ij[3],ij[4],
                                     json$rois[[ij[1]]],
                                     json$tracks[[ij[1]]])
                    ij[1:2] <- c(ij[1]+1,1)
                    ij <- load.grain(dl,ij,mica,json,new=FALSE,lab=lab)
                } else if (donew(xy)){
                    #nrois <- length(json$rois)
                    #jsron$rois[[nrois+1]] <- list(x=NULL,y=NULL)
                } else if (doexit(xy)){
                    if (!identical(response,'b'))
                        save.to.json(dl[ij[1]],ij[3],ij[4],
                                     json$rois[[ij[1]]],
                                     json$tracks[[ij[1]]])
                    break;
                } else if (dodelete(xy) & identical(response,'r')){
                    roi <- json$rois[[ij[1]]]
                    nc <- length(roi$x)
                    if (nc==1){
                        roi$x <- NULL
                        roi$y <- NULL
                    } else {
                        roi$x <- roi$x[1:(nc-1)]
                        roi$y <- roi$y[1:(nc-1)]
                    }
                    json$rois[[ij[1]]] <- roi
                    load.grain(dl,ij,mica,json,new=TRUE,lab=lab)
                } else if (dodelete(xy) & identical(response,'c')){
                    tracks <- json$tracks[[ij[1]]]
                    nt <- length(tracks$x)
                    if (nt==1){
                        tracks$x <- NULL
                        tracks$y <- NULL
                    } else {
                        tracks$x <- tracks$x[1:(nt-1)]
                        tracks$y <- tracks$y[1:(nt-1)]
                    }
                    json$tracks[[ij[1]]] <- tracks
                    load.grain(dl,ij,mica,json,new=TRUE,lab=lab)
                } else if (domica(xy)){
                    mica <- TRUE
                    ij[2] <- 1
                    load.grain(dl,ij,mica,json,new=FALSE,lab=lab)
                } else if (dograin(xy)){
                    mica <- FALSE
                    ij[2] <- 1
                    load.grain(dl,ij,mica,json,new=FALSE,lab=lab)
                } else if (identical(response,'r')) {
                    json$rois[[ij[1]]]$x <- c(json$rois[[ij[1]]]$x,xy$x)
                    json$rois[[ij[1]]]$y <- c(json$rois[[ij[1]]]$y,1-xy$y)
                    load.grain(dl,ij,mica,json,new=TRUE,lab=lab)
                } else if (identical(response,'c')) {
                    json$tracks[[ij[1]]]$x <- c(json$tracks[[ij[1]]]$x,xy$x)
                    json$tracks[[ij[1]]]$y <- c(json$tracks[[ij[1]]]$y,1-xy$y)
                    load.grain(dl,ij,mica,json,new=TRUE,lab=lab)
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

draw.spots <- function(tracks){
    nt <- length(tracks$x)
    if (nt>0){
        for (i in 2:nt){
            points(tracks$x[(i-1):i],1-tracks$y[(i-1):i],pch=21,col='yellow')
        }
    }
}

save.to.json <- function(dir,height,width,roi,tracks){
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
    out <- paste0(out,'],"shift":[0,0]}]')
    nt <- length(tracks$x) # number of tracks
    out <- paste0(out,',"tracks":[')
    for (j in 1:nt){
        out <- paste0(out,'[',round(width*tracks$x[j]),
                      ',',round(height*tracks$y[j]),']')
        if (j<nt) out <- paste0(out,",")    
    }
    out <- paste0(out,']}')
    sink(paste0(dir,'/rois.json'))
    cat(out)
    sink()
    if (TRUE){ # set to TRUE to feed track ROIS to tutorial
        print.json(roi,tracks,width,height,12,12)
    }
}

print.json <- function(roi,tracks,w,h,dx,dy){
    out <- '[['
    if (!is.null(roi$x)){
        for (i in 1:length(roi$x)){
            x <- round(w*roi$x[i])
            out <- paste0(out,x,',')
        }
        out <- paste0(out,round(w*roi$x[1]))
    }
    out <- paste0(out,'],[')
    if (!is.null(roi$y)){
        for (i in 1:length(roi$y)){
            y <- round(h*roi$y[i])
            out <- paste0(out,y,',')
        }
        out <- paste0(out,round(h*roi$y[1]))
    }
    out <- paste0(out,']],[[')
    for (i in 1:length(tracks$x)){
        if (!is.null(tracks$x) & !is.null(tracks$y)){
            tx <- round(w*tracks$x[i])
            ty <- round(h*tracks$y[i])
            out <- paste0(out,'[[',tx-dx,',',tx-dx,',',tx+dx,',',tx+dx,',',tx-dx)
            out <- paste0(out,'],[',ty-dy,',',ty+dy,',',ty+dy,',',ty-dy,',',ty-dy,']],')
        }
    }
    out <- paste0(out,']')
    print(out)
}

prepareGentFiles <- function(dir,ngrains){
    files <- list.files(dir)
    ndig <- nchar(ngrains)
    for (i in 1:ngrains){
        gnum <- paste0(rep(0,ndig-nchar(i)),i)
        gname <- paste0('xy',gnum)
        sfiles <- files[grepl(pattern=gname,files)]
        system(paste0('mkdir ',dir,'Grain',gnum))
        cat(paste0('{"image_width":1608,"image_height":1608,',
                   '"regions":[{"vertices":[],"shift":[]}]}'),
            file=paste0(dir,'Grain',gnum,'/rois.json'))
        for (file in sfiles){
            system(paste0('mv ',dir,file,' ',dir,'Grain',gnum))
        }
    }
}

#prepareGentFiles('/home/pvermees/Desktop/GaH/FCT-G1/',62)

#geochron.at.home("LU324-2-DUR")
#geochron.at.home("LU324-6-FCT")
#geochron.at.home("/home/pvermees/Documents/fissiontracks/geochron@home/images/LU324-2-DUR")
#geochron.at.home("/home/pvermees/Documents/geochron@home/LU324-10-MD")
#geochron.at.home("/home/pvermees/Documents/geochron@home/LU288-2-GIL")
#geochron.at.home("/home/pvermees/Documents/fissiontracks/geochron@home/images/Qingyang-1X")

geochron.at.home("/home/pvermees/Desktop/GaH/FCT-G1",lab='Gent')
