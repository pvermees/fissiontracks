rm(list=ls())

xml2json <- function(fname,odir,sample="mysample",project="myproject"){
    
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
        rois[[i]] <- parseROI(f)
        results[[i]] <- grain2grain(g,sample,project)
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
        raw <- jsonlite::toJSON(results[[i]],auto_unbox=TRUE)
        opening <- gsub('"results":{',replacement='"results":[{',x=raw,fixed=TRUE)
        json <- gsub(']]',replacement=']}],"area_pixels"',x=opening,fixed=TRUE)
        cat(json,file=file.path(oname,"results.json"))
    }
}

parseROI <- function(f){
    list(
        image_width = as.numeric(f$intGrainWidth),
        image_height = as.numeric(f$intGrainHeight),
        regions = parseVertices(f$roiRegionOfInterest)
    )
}

parseVertices <- function(roi){
    nparts <- lengths(regmatches(roi, gregexpr("#", roi)))
    parts <- strsplit(roi,split="#")[[1]]
    out <- list()
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
    
grain2grain  <- function(g,sample,project){
    stagexyz <- strsplit(g$Fields$locGrainLocation,",")[[1]]
    list(
        grain = 0,
        sample = 0,
        sample_name = sample,
        project_name = project,
        results = list(
            ft_type = "S",
            result = length(g$Tracks),
            create_date = date(),
            worker = list(id=0),
            latlngs = parseTracks(g$Tracks),
            area_pixels = 0798323.5,
            area_mm2 = 0.006581889861039999,
            id = 0,
            index = 0,
            image_width = as.numeric(g$Fields$intGrainWidth),
            image_height = as.numeric(g$Fields$intGrainHeight),
            scale_x = 9.079999999999999e-08,
            scale_y = 9.079999999999999e-08,
            stage_x = as.numeric(stagexyz[2]),
            stage_y = as.numeric(stagexyz[3]),
            mica_stage_x = NULL,
            mica_stage_y = NULL,
            shift_x = 0,
            shift_y = 0
        )
    )
}

parseTracks <- function(tracks){
    nt <- length(tracks)
    mat <- matrix(0,nt,2)
    for (i in seq_along(tracks)){
        mat[i,] <- as.numeric(strsplit(tracks[[i]]$Fields$pntCenter,",")[[1]])
    }
    out <- asplit(mat,1)
    out
}

tif2jpeg <- function(idir,odir){
    grains <- list.files(idir,pattern="Grain*")
    for (grain in grains){
        message(grain)
        igrain <- file.path(idir,grain)
        ograin <- file.path(odir,grain)
        if (!dir.exists(ograin)) dir.create(ograin)
        itrans <- file.path(igrain,"Stack.tif")
        otrans <- file.path(ograin,"Stack.jpg")
        system(paste0("convert -set colorspace Gray -quality 50 ",itrans," ",otrans))
        iflat <- file.path(igrain,"ReflStackFlat.tif")
        oflat <- file.path(ograin,"ReflStackFlat.jpg")
        system(paste0("convert -set colorspace Gray  -quality 50 ",iflat," ",oflat))
    }
}

FT2GaH <- function(idir,odir,xml){
    xml2json(fname=file.path(idir,xml),sample="mysample",
             project="myproject",odir=odir)
    tif2jpeg(idir=idir,odir=odir)
}

if (TRUE) { # example:
    FT2GaH(idir="~/Documents/FTrawData/20HL-06-FT",
           odir="~/Documents/FTrawData/20HL-06-GaH",xml="20HL-06.xml")
}
