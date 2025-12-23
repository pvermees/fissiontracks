options(warn=0)

get_json <- function(image_width=2736.0,
                     image_height=2208.0,
                     x,y){
    paste0('{"image_width": ',image_width,
           ', "image_height": ',image_height,
           ', "scale_x": 9.08e-08, "scale_y": 9.08e-08, ',
           '"stage_x": ',x,', "stage_y": ',y,
           ', "regions": [{"vertices": [], "shift": [0, 0]}]}')
}

mixedsort <- function(grains){
    extracted_numbers_char <- gsub(".*\\((\\d+)\\).*", "\\1", grains)
    extracted_numbers <- as.numeric(extracted_numbers_char)
    sort_index <- order(extracted_numbers, na.last = FALSE)
    grains[sort_index]
}

MiA2GaH <- function(idir,odir,IMpath=NULL){
    grains <- list.files(idir,pattern="Grain*")
    convert <- ifelse(is.null(IMpath),"convert",file.path(IMpath,"convert"))
    poi <- read.csv('poi.csv',check.names=FALSE)
    x <- poi[,"x::x"]
    y <- poi[,"y::y"]
    sortedgrains <- mixedsort(grains)
    for (i in seq_along(sortedgrains)){
        print(i)
        grain <- grains[i]
        igrain <- file.path(idir,grain)
        images <- list.files(igrain,pattern="\\.png$")
        if (i<10) gname <- paste0("Grain0",i)
        else gname <- paste0("Grain",i)
        gdir <- file.path(odir,gname)
        if (!dir.exists(gdir)) dir.create(gdir)
        for (j in seq_along(images)){
            ifname <- file.path(idir,paste0("'",grain,"'"),images[j])
            if (j<10) stackname <- paste0("Stack-0",j)
            else stackname <- paste0("Stack-",j)
            ofname <- file.path(gdir,paste0(stackname,".jpg"))
            cmd <- paste0(convert," -set colorspace Gray -quality 50 ",ifname," ",ofname)
            system(cmd)
        }
        cat(get_json(x=x[i],y=y[i]),file=file.path(gdir,'rois.json'))
    }
}

setwd("~/Documents/fissiontracks/geochron@home/Navya_Veepuru/DEN71-Ind1835")

MiA2GaH(idir="Outputs",odir="GaH")
