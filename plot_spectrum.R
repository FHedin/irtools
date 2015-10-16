#!/usr/bin/env Rscript

# args <- commandArgs(trailingOnly = TRUE)

# inp = read.table(args[1])
inp = read.table("spectrum.cpu")

mat = matrix(inp$V3,nrow=sqrt(length(inp$V1)),ncol=sqrt(length(inp$V2)),byrow=FALSE) 

themin = min(mat) 
themax = max(mat)

mat=mat/themax

themin = min(mat) 
themax = max(mat)

vx=c(1000,1065)
vy=c(1000,1065)

ym = seq(min(inp$V1),max(inp$V1),length.out=nrow(mat))
xm = seq(min(inp$V2),max(inp$V2),length.out=ncol(mat))

# png('spec2d.png',width=4000,height=4000,pointsize=96)

# par(bty="o",pin=c(7.0,7.0),mai=c(8.2,8.2,1.2,1.2))

# veclmi = seq(from=themin,to=-0.1,by=0.1)
# veclma = seq(from=0.1,to=themax,by=0.1)
veclmi = seq(from=-1.0,to=-0.1,by=0.1)
veclma = seq(from=0.1,to=1.0,by=0.1)


xlb = expression(omega[3]*" "*(cm^-1))
ylb = expression(omega[1]*" "*(cm^-1))

contour(xm,ym,mat,levels=veclmi,asp=1.0,col="red",xlab=xlb,ylab=ylb,drawlabels=T)  #,lwd=9,cex.lab=2,cex.axis=2)

contour(xm,ym,mat,levels=veclma,asp=1.0,col="blue",add=TRUE,drawlabels=T)    #,lwd=9,cex.lab=2,cex.axis=2)


# axis(1,lwd=9,cex.lab=2,cex.axis=2,lwd.ticks=9)
# axis(2,lwd=9,cex.lab=2,cex.axis=2,lwd.ticks=9)
# 
# lines(xm,ym,col="black",lwd=9)
# 
# box(lwd=9)

# dev.off()

