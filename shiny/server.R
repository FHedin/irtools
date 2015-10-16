library(shiny)

shinyServer(function(input, output) {
  
  output$plot2dir <- renderPlot({
    
    file = gzfile("test_spectrum.dat.gz")
    inp = read.table(file)
    
    mat = matrix(inp$V3,nrow=sqrt(length(inp$V1)),ncol=sqrt(length(inp$V2)),byrow=FALSE) 
    
    themin = min(mat) 
    themax = max(mat)
    
    mat=mat/themax
    
    themin = min(mat) 
    themax = max(mat)
    
    vx=input$axRange 
    vy=input$axRange
    
    ym = seq(min(inp$V1),max(inp$V1),length.out=nrow(mat))
    xm = seq(min(inp$V2),max(inp$V2),length.out=ncol(mat))

    veclmi = seq(from=-1.0,to=-0.1,by=0.1)
    veclma = seq(from=0.1,to=1.0,by=0.1)
    
    xlb = expression(omega[3]*" "*(cm^-1))
    ylb = expression(omega[1]*" "*(cm^-1))
    par(pty="s",lwd=2,cex=1.5)
    contour(xm,ym,mat,levels=veclmi,
            col="red",xlab=xlb,ylab=ylb,drawlabels=T,
            xlim=vx, ylim=vy,asp=1,axes=F)
    contour(xm,ym,mat,levels=veclma,col="blue",add=TRUE,drawlabels=T)   
    
     axis(1)
     axis(2)
    
    lines(xm,ym,col="black")
    
    box()
    
  })
  
})
