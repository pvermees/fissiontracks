options(warn=2)

MiA2GaH <- function(idir,odir,IMpath=NULL){
    grains <- list.files(idir,pattern="Grain*")
    convert <- ifelse(is.null(IMpath),"convert",file.path(IMpath,"convert"))
    for (grain in grains){
        message(grain)
        igrain <- paste0(idir,grain)
        ograin <- paste0(odir,grain)
        if (!dir.exists(ograin)) dir.create(ograin)
        images <- list.files(igrain,pattern="\\.png$")
        for (i in seq_along(images)){
            ifname <- file.path(igrain,images[i])
            ofname <- file.path(ograin,paste0("Stack-0",i-1,".jpg"))
            cmd <- paste0(convert," -set colorspace Gray -quality 50 '",ifname,"' '",ofname,"'")
            system(cmd)
        }
        foo <- 1
    }
}

MiA2GaH(idir="~/Documents/fissiontracks/geochron@home/AndyCarter/SEY809/Outputs/",
        odir="~/Documents/fissiontracks/geochron@home/AndyCarter/SEY809/GaH/")
