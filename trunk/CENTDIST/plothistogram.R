#to run:
#R input w binsize output --no-save --args --slave < plothistogram.R
Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]
plothistogram=function(x,left,right,binsize,oddbins=T,col="blue",xlab="xlab",ylab="ylab",main="main"){
center=(left+right)%/%2
if (oddbins){
w=((right-left)%/%2-binsize/2)%/%binsize*binsize+binsize/2
hist(x[-w<=x-center & x-center<w],breaks=(((center-w-binsize/2)/binsize):((center+w-binsize/2)%/%binsize))*binsize+binsize/2,col=col,xlab=xlab,ylab=ylab,main=main)
}
else{
w=((right-left)%/%2)%/%binsize*binsize
hist(x[-w<=x-center & x-center<w],breaks=(((center-w)/binsize):((center+w)%/%binsize))*binsize,col=col,xlab=xlab,ylab=ylab,main=main)
}
}
x<-read.table(Args[2])
w<-Args[3]
binsize<-Args[4]
jpeg(Args[5])
plothistogram(x[,2],-w,w,binsize)
dev.off()
