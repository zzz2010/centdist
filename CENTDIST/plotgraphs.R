plothistogram=function(x,left,right,binsize,oddbins=T,col="blue",xlab="xlab",ylab="ylab",main="main"){
center=(left+right)%/%2
if (oddbins){
w=((right-left)%/%2-binsize/2)%/%binsize*binsize+binsize/2
hist(x[-w<=x-center & x-center<w],breaks=(((center-w-binsize/2)%/%binsize):((center+w-binsize/2)%/%binsize))*binsize+binsize/2,col=col,xlab=xlab,ylab=ylab,main=main)
}   
else{
w=((right-left)%/%2)%/%binsize*binsize
hist(x[-w<=x-center & x-center<w],breaks=(((center-w)%/%binsize):((center+w)%/%binsize))*binsize,col=col,xlab=xlab,ylab=ylab,main=main)
}
}
plotfoldhistogram=function(x,left,right,binsize,col="blue",xlab="xlab",ylab="ylab",main="main"){
center=(left+right)%/%2
w=((right-left)%/%2)%/%binsize*binsize
hist( abs(x[-w<=x-center & x-center<w]),breaks=(0:((center+w)%/%binsize))*binsize,col=col,xlab=xlab,ylab=ylab,main=main)
}

Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]
inputdir<-Args[2] 
readdict=function(f){ x<-read.table(f,as.is=T) 
	y<-x[,2]
	names(y)<-x[,1]
	y
}
wd<-getwd()
setwd(inputdir)
x<-read.table('info.txt',as.is=T)
info<-x[,2]
names(info)<-x[,1]
N<-as.integer(info['N'])
w<-as.integer(info['w'])
L<-as.integer(info['L'])
read.table('results.txt',header=T)->results
results<-results[order(results[,"SCORE"],decreasing=T),]
matlens<-readdict('matlist.txt')

for (i in 1:length(results[,1])){
print(i)
motifname<-as.character(results[i,"NAME"])
motifname<-basename(substring(motifname,1,nchar(motifname)-4))
posfile<-paste(motifname,".pos",sep="")
x<-read.table(posfile)
x[,2]=x[,2]+matlens[motifname]%/%2-L%/%2
y<-x[x[,3]>=results[i,"THRESHOLD"],]
graphics.off()
#jpeg(paste(posfile,"disthist.jpg",sep="_"));
#plothistogram(y[,2],-L/2,L/2,50,xlab="",ylab="",main=motifname)
#dev.off()
#graphics.off()
jpeg(paste(posfile,"disthist.jpg",sep="_"),width=200,height=200)
par(mar=c(2,2,2,1),cex.axis=0.9)
plothistogram(y[,2],-L/2,L/2,50,xlab="",ylab="",main=motifname)
dev.off()
graphics.off()
#jpeg(paste(posfile,"distfoldhist.jpg",sep="_"));
#plotfoldhistogram(y[,2],-L/2,L/2,50,xlab="",ylab="",main=motifname)
#dev.off()
#graphics.off()
jpeg(paste(posfile,"distfoldhist.jpg",sep="_"),width=200,height=200)
par(mar=c(2,2,2,1),cex.axis=0.9)
plotfoldhistogram(y[,2],-L/2,L/2,50,xlab="",ylab="",main=motifname)
dev.off()

z=y[abs(y[,2])<results[i,"W"],1]
z=c(0,z[!duplicated(z)]+1,N+1)
numhits=length(z)
xindex=1:numhits
xaxis=z[xindex];
n=pmin(numhits%/%4,100)
right=pmin(numhits,xindex+n)
left=pmax(1,xindex-n)
yaxis=(right-left)/(z[right]-z[left])
if (n>=5){
jpeg(paste(posfile,"movingaverage.jpg",sep="."),height=200,width=200);
par(mar=c(2.5,2,2,0.4),cex.axis=0.9,cex.lab=0.9)
plot(xaxis,yaxis,xlab='Peak Rank',ylab='% Sequence containing motif',main='',ylim=c(0,max(yaxis)),type='l',col='blue',mgp=c(1.2,0.5,0))
dev.off()
}
z=x[abs(x[,2])<200,1]
z=c(0,z[!duplicated(z)]+1,N+1)
numhits=length(z)
xindex=1:numhits
xaxis=z[xindex];
n=pmin(numhits%/%4,100)
right=pmin(numhits,xindex+n)
left=pmax(1,xindex-n)
yaxis=(right-left)/(z[right]-z[left])
if (n>=5){
jpeg(paste(posfile,"movingaverage2.jpg",sep="."),height=250,width=250);
par(mar=c(2.5,2,2,0.4),cex.axis=0.9,cex.lab=0.9)
plot(xaxis,yaxis,xlab='Peak Rank',ylab='% Sequence containing motif',main='',ylim=c(0,max(yaxis)),type='l',col='blue',mgp=c(1.2,0.5,0))
dev.off()
}
Z=c(-w,y[,2],w);
tableZ=as.matrix(table(Z));
z=1:sum(tableZ)
k=0
for (i2 in row.names(tableZ)){
	count=tableZ[i2,]
	if (count>0){
		for (j in 1:count){
			k=k+1
			z[k]=as.numeric(i2)+j/(count+1)
		}
	}
}
numhits=length(z)
xindex=1:numhits
xaxis=z[xindex];
n=pmin(numhits%/%4,100)
right=pmin(numhits,xindex+n)
left=pmax(1,xindex-n)
yaxis=(right-left)/(z[right]-z[left])
if (n>=5){
jpeg(paste(posfile,"disthist.new.jpg",sep="."),width=250,height=250);
par(mar=c(2.5,2,2,0.4),cex.axis=0.9,cex.lab=0.9)
plot(xaxis,yaxis,xlab='',ylab='',main='',ylim=c(0,max(yaxis)),type='l',col='blue',mgp=c(1.2,0.5,0))
dev.off()
}
z=y[abs(y[,2])<results[i,"W"],1]
z=c(z[!duplicated(z)]+1)
numhits=length(z)
write(numhits/N,paste(motifname,".count1",sep=""))
z=x[abs(x[,2])<200,1]
z=c(z[!duplicated(z)]+1)
numhits=length(z)
write(numhits/N,paste(motifname,".count2",sep=""))
#break;
}
