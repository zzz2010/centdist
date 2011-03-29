#to run:
#R input w binsize output --no-save --args --slave < plotmovingaverage.R
Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]
f=Args[2]
N=as.numeric(Args[3])
binsize=as.numeric(Args[4])
out=paste(f,'.movingaverage.jpg',sep='')
movingaverage=function(x,n){
	cx=c(0,cumsum(x))
	N=length(x)
	(cx[2:(N+1)]-cx[pmax(2:(N+1)-n,1)])/pmin(1:N,n)
	#c(rep(NA,n-1),cx[n:N+1]-cx[n:N+1-n])/n
}
plotmovingaverage=function(x,N,binsize){
	#for plotting moving average of the motif occurrence
	m=movingaverage(table(c(0:(N-1),x))-table(0:(N-1))>0,n=binsize)
	plot(1:N,m,ylab="Percentage sequence containing motif",xlab="Peak Rank",main='',ylim=c(0,max(m,na.rm=T)),type='l',col='blue')
}
x=read.table(f);
L=max(x[,2])
y=x[abs(x[,2]-L/2)<200,1]
y=c(0,y[!duplicated(y)]+1)
N=length(y)
n=pmin(N%/%2,100)
xaxis=y[(((n+2):N)+((n+2):N-(n+1)))%/%2]
yaxis=n/(y[(n+2):N]-y[(n+2):N-(n+1)])
graphics.off();jpeg(out);
#plotmovingaverage(x[abs(x[,2]-L/2)<200,1],N,binsize);
plot(xaxis,yaxis,xlab='Peak Rank',ylab='Percentage Sequence containing motif',main='',ylim=c(0,max(yaxis)),type='l',col='blue')
dev.off()
