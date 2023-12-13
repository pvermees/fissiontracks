xml2json <- function(fname,odir,user="admin",sample="mysample"){
    
    rois <- list()
    counts <- list()

    txt <- readChar(fname,file.info(fname)$size)
    xml <- XML::xmlParse(txt)
    message("Converting the XML file to a list.")
    XMLlist <- XML::xmlToList(xml)
    message("Converting the list to JSON.")

    ng <- length(XMLlist$Grains)
    XMLgrains <- XMLlist$Grains

    rois <- list()
    results <- list()
    for (i in seq_along(XMLgrains)){
        g <- XMLgrains[[i]]
        f <- g$Fields
        rois[[i]] <- parseROI(f,debug=(i==12))
        if (is.null(rois[[i]]$regions)){
            warning('grain ',i,' does not have a region of interest')
        }
        results[[i]] <- grain2grain(g,user=user,sample=sample,index=i)
    }

    for (i in seq_along(rois)){
        padding <- ifelse(i<10,"0","")
        dname <- paste0("Grain",padding,i)
        oname <- file.path(odir,dname)
        if (!dir.exists(oname)) dir.create(oname)
        raw <- jsonlite::toJSON(rois[[i]],auto_unbox=TRUE)
        opening <- gsub('"regions":{',replacement='"regions":[{',x=raw,fixed=TRUE)
        json <- gsub('}}',replacement='}]}',x=opening,fixed=TRUE)
        cat(json,file=file.path(oname,"rois.json"))
    }
    json <- jsonlite::toJSON(results,auto_unbox=TRUE)
    cat(json,file=file.path(odir,"results.json"))

}

parseROI <- function(f,debug=FALSE){
    list(
        image_width = as.numeric(f$intGrainWidth),
        image_height = as.numeric(f$intGrainHeight),
        regions = parseVertices(f$roiRegionOfInterest)
    )
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
    
grain2grain  <- function(g,user,sample,index){
    stagexyz <- strsplit(g$Fields$locGrainLocation,",")[[1]]
    list(
        sample = sample,
        index = index,
        user = user,
        date = date(),
        ft_type = "S",
        points = parseTracks(g$Tracks)
    )
}

parseTracks <- function(tracks){
    nt <- length(tracks)
    out <- list()
    for (i in seq_along(tracks)){
        xy <- as.numeric(strsplit(tracks[[i]]$Fields$pntCenter,",")[[1]])
        out[[i]] <- list(x_pixels=xy[1],y_pixels=xy[2])
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
        oflat <- file.path(ograin,"ReflStackFlat.jpg")
        system(paste0(convert," -set colorspace Gray  -quality 50 ",iflat," ",oflat))
    }
}

FT2GaH <- function(idir,odir,xml,user="admin",sample,IMpath=NULL){
    xml2json(fname=file.path(idir,xml),odir=odir,user=user,sample=sample)
    tif2jpeg(idir=idir,odir=odir,IMpath=IMpath)
}
