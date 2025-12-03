options(warn=2)

MiA2GaH <- function(idir,odir,IMpath=NULL){
    grains <- list.files(idir,pattern="Grain*")
    convert <- ifelse(is.null(IMpath),"convert",file.path(IMpath,"convert"))
    for (i in seq_along(grains)){
        grain <- grains[i]
        message(grain)
        igrain <- paste0(idir,grain)
        images <- list.files(igrain,pattern="\\.png$")
        if (i<10) gname <- paste0("Grain0",i)
        else gname <- paste0("Grain",i)
        gdir <- paste0(odir,gname)
        if (!dir.exists(gdir)) dir.create(gdir)
        for (j in seq_along(images)){
            ifname <- file.path(paste0(idir,"'",grain,"'"),images[j])
            if (j<10) stackname <- paste0("Stack-0",j)
            else stackname <- paste0("Stack-",j)
            ofname <- file.path(gdir,paste0(stackname,".jpg"))
            cmd1 <- paste0(convert," -set colorspace Gray -quality 50 ",ifname," ",ofname)
            ixmlname <- file.path(paste0(idir,"'",grain,"'"),paste0(images[j],"_metadata.xml"))
            oxmlname <- file.path(gdir,paste0(stackname,".jpg_metadata.xml"))
            cmd2 <- paste0("cp ",ixmlname," ",oxmlname)
            system(cmd1)
            system(cmd2)
        }
        foo <- 1
    }
}

MiA2GaH(idir="~/Documents/fissiontracks/geochron@home/AndyCarter/SEY809/Outputs/",
        odir="~/Documents/fissiontracks/geochron@home/AndyCarter/SEY809/GaH/")
