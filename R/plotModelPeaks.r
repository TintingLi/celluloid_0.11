

# epp is obtained from a call to ePeakPos 
# ... passed to addLabels 
plotModelPeaks<-function( par=NULL, S=NULL ,t=NULL , selectedPoints=NULL , cn, epcol="red", epcex=3, eplwd=3 , eppch=21, 
                          spcol="black", sppch=19, spcex=1, splwd=1, addlabels=F,  preserveMatPatDiff=T , preserveMaxSubClDiff=T , 
                          plot=TRUE, xonly=FALSE, ... ){
  
  if( !is.null(par) ){
    if( !is.null( S ) | !is.null(t) ){
      stop("only need one of par or S/t")
    }
    S<-par[1]
    t<-par[2:length(par)]
    t<-c(t, 1-sum(t))
  } else {
    if( is.null(S) | is.null(t) ){
      stop("need either par or S/t")
    }
  }
  
  if( !xonly ){
    epp<- ePeakPos( S=S, t=t, cn=cn, preserveMatPatDiff=preserveMatPatDiff, preserveMaxSubClDiff=preserveMaxSubClDiff   )
  } else {
    epp<- ePeakPos.xonly( S=S, t=t, cn=cn ) 
  }
  
  
  if( plot ){
    points( unique( epp[,c("x","ar")] ),  pch=eppch , col=epcol , cex=epcex, lwd=eplwd )
    
    # draw vertical line where all tumour cells have the same integer copy number count 
    sel<- apply( epp, 1, function(x){  return( all( x[seq( 3,length(x)-2, 2 ) ]==x[3] &  x[seq( 4,length(x)-2, 2 ) ]==x[4]  ) ) } ) 
    
    integercn <- unique( apply( epp[sel,3:4], 1, sum) )
    x<-unique( epp[sel,"x"] )
    abline( v=x, col=epcol)
    #  axis( 1, at=x, labels= sort( unique( integercn ) ) , col.axis=epcol ) 
    
    #this won't work if the subclones are in equal freq
    #ep<- unique( apply( cn, 1, RCC, S=S, t=t ) ) 
    ep<- ( apply( cn, 1, celluloid:::RCC, S=S, t=t ) ) 
    
    arp<-prepAR( t )
    d<-data.frame(cn,ep)
    d<-d[ order( d$ep),] 
    
    s<-seq(1,length(ep),2)
    axis( 1, at=d$ep[s] , labels=apply( as.matrix(d[,2:ncol(cn)]), 1, paste, collapse =".")[s] , las=2  , cex.axis=.75) 
    s<-seq(2,length(ep),2)
    axis( 3, at=d$ep[s] , labels=apply( as.matrix(d[,2:ncol(cn)]), 1, paste, collapse =".")[s] , las=2  , cex.axis=.75  ) 
    
    
    if( !is.null(selectedPoints) ){
      points( selectedPoints[,1:2], pch=sppch,  col=spcol , cex=spcex, lwd=splwd )
      
      if( addlabels)
        addLabels( epp, selectedPoints, ... )
      
    }
    
    param<-floor( 1000*t+.5 )/1000
    pl<-(2/S-2*param[1])/(1-param[1])
    
    tit<-paste( c(paste("S:",floor( 1000*S+.5)/1000, sep=""), 
                  paste("PLOIDY:", floor( 1000*pl )/1000 , sep=""),
                  paste( "%N:",param[1], sep=""),
                  paste( paste( "%T",1:(length(param)-1) , sep=""  ), param[2:length(param)], sep=":" )  ), collapse=", ")
    
    title( main=tit , cex.main= 1.25 - length(param)/12 , line= 2.5   )
  }
  invisible(epp) 
  
}



