xml2json <- function(fname,odir,user,analyst,sample){
    
    rois <- list()
    counts <- list()

    txt <- readChar(fname,file.info(fname)$size)
    xml <- XML::xmlParse(txt)
    message("Converting the XML file to a list.")
    XMLlist <- XML::xmlToList(xml)
    message("Converting the list to JSON.")

    scaleZ <- as.numeric(XMLlist$AllData$dblTransStepSize)

    ng <- length(XMLlist$Grains)
    XMLgrains <- XMLlist$Grains

    rois <- list()
    results <- list()
    for (i in seq_along(XMLgrains)){
        # 1. track counts and lengths
        g <- XMLgrains[[i]]
        f <- g$Fields
        gname <- f$strGrainName
        rois[[i]] <- parseROI(f,debug=(i==12))
        gnum <- as.numeric(gsub("\\D", "", gname))
        results[[i]] <- grain2grain(g,user=user,analyst=analyst,
                                    sample=sample,index=gnum,scaleZ=scaleZ)
        # 2. rois
        regions <- rois[[i]]$regions
        if (!is.null(names(regions))){ # turn vector into list
            rois[[i]]$regions <- list()
            rois[[i]]$regions[[1]] <- regions
        }
        results[[i]]$roi <- rois[[i]]
        gname <- XMLgrains[[i]]$Fields$strGrainName
        oname <- file.path(odir,gname)
        if (!dir.exists(oname)) dir.create(oname)
        json <- jsonlite::toJSON(rois[[i]],auto_unbox=TRUE)
        cat(json,file=file.path(oname,"rois.json"))
    }

    json <- jsonlite::toJSON(results,auto_unbox=TRUE)
    cat(json,file=file.path(odir,"results.json"))

}

parseROI <- function(f,debug=FALSE){
    out <- list(
        image_width = as.numeric(f$intGrainWidth),
        image_height = as.numeric(f$intGrainHeight),
        regions = parseVertices(f$roiRegionOfInterest)
    )
    if (is.null(out$regions)){
        out$regions <- list(vertices=list(c(0,0),
                                          c(out$image_width,0),
                                          c(out$image_width,out$image_height),
                                          c(0,out$image_height)),
                            shift=c(0,0))
    }
    out
}

parseVertices <- function(roi){
    if (is.null(roi)) return(NULL)
    out <- list()
    nparts <- lengths(regmatches(roi, gregexpr("#", roi)))
    parts <- strsplit(roi,split="#")[[1]]
    for (i in 1:nparts){
        part <- gsub("true;","",parts[i])
        v1 <- unlist(strsplit(part,split=";"))
        v2 <- as.numeric(unlist(strsplit(v1,split=",")))
        nr <- length(v2)
        mat <- cbind(v2[seq(from=1,to=nr-1,by=2)],
                     v2[seq(from=2,to=nr,by=2)])
        out[[i]] <- list(vertices=asplit(mat,1), shift=c(0,0))
    }
    out
}
    
grain2grain  <- function(g,user,analyst,sample,index,scaleZ){
    stagexyz <- strsplit(g$Fields$locGrainLocation,",")[[1]]
    list(
        sample = sample,
        index = index,
        user = user,
        analyst = analyst,
        date = date(),
        ft_type = "S",
        points = parseTracks(g$Tracks),
        lines = parseLengths(g$Fields$tvlTINTs3d,scaleZ=scaleZ)
    )
}

parseTracks <- function(tracks){
    out <- list()
    if (is.null(tracks)) return(out)
    else nt <- length(tracks)
    j <- 0
    for (i in seq_along(tracks)){
        ai <- "Outline" %in% names(tracks[[i]])
        fields <- tracks[[i]]$Fields
        if (isa(fields,'list')){
            nt <- fields$intMultiple
            xy <- as.numeric(strsplit(fields$pntCenter,",")[[1]])
        } else {
            fieldsvec <- strsplit(fields,";")[[1]]
            if (ai){
                nt <- as.numeric(fieldsvec[4])
                xy <- as.numeric(strsplit(fieldsvec[10],"[:,]")[[1]][2:3])
            } else {
                nt <- as.numeric(fieldsvec[3])
                xy <- as.numeric(strsplit(fieldsvec[1],",")[[1]])
            }
        }
        j <- j+1
        out[[j]] <- list(x_pixels=xy[1],y_pixels=xy[2])
        if (nt>1){
            j <- j+1
            out[[j]] <- list(x_pixels=xy[1]+1,y_pixels=xy[2]+1)
        }
        if (nt>2){
            j <- j+1
            out[[j]] <- list(x_pixels=xy[1]+1,y_pixels=xy[2]-1)
        }
        if (nt>3){
            j <- j+1
            out[[j]] <- list(x_pixels=xy[1]-1,y_pixels=xy[2]-1)
        }
        if (nt>4){
            j <- j+1
            out[[j]] <- list(x_pixels=xy[1]-1,y_pixels=xy[2]+1)
        }
        if (nt>5){
            for (ii in 6:nt){
                j <- j+1
                xyj <- jitter(xy)
                out[[j]] <- list(x_pixels=xyj[1],y_pixels=xyj[2])
            }
        }
    }
    out
}

parseLengths <- function(lengths,scaleZ=1){
    out <- list()
    if (!is.null(lengths)){
        string <- substr(lengths,1,nchar(lengths)-2)
        blocks <- strsplit(string,':')[[1]]
        for (i in seq_along(blocks)){
            midpoints <- strsplit(blocks[i],split=';')[[1]]
            from <- as.numeric(strsplit(midpoints[1],split=',')[[1]])/c(1,1,scaleZ)
            to <- as.numeric(strsplit(midpoints[2],split=',')[[1]])/c(1,1,scaleZ)
            out[[i]] <- round(c(from,to))
        }
    }
    out
}

tif2jpeg <- function(idir,odir,IMpath=NULL){
    grains <- list.files(idir,pattern="Grain*")
    for (grain in grains){
        message(grain)
        igrain <- file.path(idir,grain)
        ograin <- file.path(odir,grain)
        if (!dir.exists(ograin)) dir.create(ograin)
        itrans <- file.path(igrain,"Stack.tif")
        otrans <- file.path(ograin,"Stack.jpg")
        convert <- ifelse(is.null(IMpath),"convert",file.path(IMpath,"convert"))
        system(paste0(convert," -set colorspace Gray -quality 50 ",itrans," ",otrans))
        iflat <- file.path(igrain,"ReflStackFlat.tif")
        if (file.exists(iflat)){
            oflat <- file.path(ograin,"ReflStackFlat.jpg")
        } else {
            iflat <- file.path(igrain,"ReflStack.tif[0]")
            oflat <- file.path(ograin,"ReflStackFlat.jpg")
        }
        system(paste0(convert," -set colorspace Gray  -quality 50 ",iflat," ",oflat))
    }
}

FT2GaH <- function(idir,odir,xml,user="admin",analyst=user,sample,IMpath=NULL,tif2jpg=TRUE){
    if (!dir.exists(odir)) dir.create(odir)
    xml2json(fname=file.path(idir,xml),odir=odir,user=user,analyst=analyst,sample=sample)
    if (tif2jpg) tif2jpeg(idir=idir,odir=odir,IMpath=IMpath)
}
