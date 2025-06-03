# microscope_coordinates is a .csv file formatted like this:
# sname,x_zeiss,y_zeiss
# 4660,49374.514,16888.863
# 4661,37741.731,6029.456
# 4662,42144.836,11402.271
# ...
# tie_points is a .csv file formatted like this:
# sname,x_icp,y_icp,z_icp
# 4678,57.658,43.705,2.502
# 4679,56.404,50.468,2.502
# 4680,49.199,43.491,2.502
# ...
coco <- function(microscope_coordinates,tie_points){
    xy <- data.matrix(read.csv(microscope_coordinates,header=TRUE,row.names=1))
    xyz <- data.matrix(read.csv(tie_points,header=TRUE,row.names=1))
    snames <- rownames(xy)
    tiepointnames <- rownames(xyz)
    from <- xy[tiepointnames,]
    to <- xyz[tiepointnames,]
    ns <- min(1000,length(snames))
    Hz <- 10
    spotsize <- 25
    headers <- c('Pattern Type','Pattern #','AutoFocus','Caption','GridSpacing','RasterSpacing',
                 'CurvePoints','ImageNo','Pattern Color','Pattern Options','CoaxPercent',
                 'PolPercent','RingPercent','TransPercent','ZoomPercent','','Pre-Ablation Pass',
                 'DepthPerPass','Dwell Time','Interstitial Pause','Energy','PassEnabled','Pass Count',
                 'Scan Speed','Rate','Spot Size','Fire Mode','Move Z','Same Direction','Burst Count',
                 'Beam Mode','XYRShutter Enabled','XYRShutter AutoRotate','XYRShutter Swap',
                 'XYRShutter Angle','XYRShutter Height','XYRShutter Width','','Ablation Pass',
                 'DepthPerPass','Dwell Time','Interstitial Pause','Energy','PassEnabled',
                 'Pass Count','Scan Speed','Rate','Spot Size','Fire Mode','Move Z',
                 'Same Direction','Burst Count','Beam Mode','XYRShutter Enabled','XYRShutter AutoRotate',
                 'XYRShutter Swap','XYRShutter Angle','XYRShutter Height','XYRShutter Width','X','Y','Z')
    template <- c(6,0,0,'Spot #1',150,55,25,0,65280,0,0,0,0,0,0,'','',
                  5,1,5,10,0,1,70,3,150,0,0,0,3,1,0,0,0,0,0,0,'','',
                  20,25,5,25,1,1,10,Hz,spotsize,0,1,0,0,1,0,0,0,0,0,0,
                  55.8773,43.0256,2.9215)
    out <- matrix(rep(template,ns),nrow=ns,byrow=TRUE)
    colnames(out) <- headers
    out[,'Pattern #'] <- 1:ns
    out[,'Caption'] <- snames[1:ns]
    tr <- provenance:::procfit(to[,-3],from)
    xyout <- tr$dilation * xy[1:ns,] %*% tr$rotation + rep(1,ns) %*% t(tr$translation)
    out[,c('X','Y')] <- xyout
    fit <- lm(Z ~ X + Y, data=data.frame(X=xyz[,1],Y=xyz[,2],Z=xyz[,3]))
    out[,'Z'] <- predict(fit, newdata=data.frame(X=xyout[,1],Y=xyout[,2]))
    route <- as.numeric(TSP::solve_TSP(TSP::TSP(dist(xyout))))
    ofname <- paste0(tools::file_path_sans_ext(microscope_coordinates),'-out.csv')
    write.csv(out[route,],file=ofname,row.names=FALSE,quote=FALSE)
}

coco('xy.csv','tie_points.csv')
