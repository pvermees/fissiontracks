rois <- IsoplotR:::fromJSON(file='roiss.json')

xy <- matrix(NA,length(rois),3)
colnames(xy) <- c('grain','x','y')

for (i in seq_along(rois)){
    grain <- rois[[i]]
    xy[i,'grain'] <- grain$grain_id
    xy[i,'x'] <- grain$stage_x
    xy[i,'y'] <- grain$stage_y
}

plot(xy[,'x'],xy[,'y'],xlab='x',ylab='y')
text(xy[,'x'],xy[,'y'],labels=xy[,'grain'],pos=1)
write.csv(xy,file='xy.csv',row.names=FALSE)
dev.copy2pdf(file='xy.pdf',width=10,height=10)
