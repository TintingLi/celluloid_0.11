
estimateLOHcurve<-function( segments , cntr=NULL, plot=F, manual=F, minsize=10000000, minmap=0.9, maxmean=1.25 ){

  ARloh<-function( x,S, n ){
    k<-2*(x-S*n)/(S-S*n) 
    return( n/( 2*n+(1-n)*k ) )
  }
  sel<-segments$size>minsize & !segments$mask & segments$meanmap>minmap  & segments$mean< maxmean & !is.na( segments$p )
  subseg <-segments[sel,]
 
  if( !manual ){ 
      f<-function( n ){
      x<-subseg$mean 
      # expected allelic ratio on the loh curve
      ear<-ARloh( x, 1, n )
      diff<- subseg$p-ear 
      diff[ diff> .025 ]<- 0 
      if( all( diff==0 ) | sum(diff!=0)<4 ){return( sum(  subseg$p-ear  )^2 ) } else { return(sum(diff^2)) }
    }
  
    Sn <- optimize( f, c(0,1) )$minimum
    
    if( plot ){
      if( is.null(cntr)  ){ stop("estimateLOHcurve: missing cntr argument")}
      image( cntr , col=terrain.colors(50) )
      contour( cntr , nlev=50, add=T ) 
      le<- subseg$end.pos-subseg$start.pos
      cxcut<- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3
      points( x<-subseg$mean, y<-subseg$p,  pch=21 , col="blue", lwd=3 , cex=cxcut  )
      points( subseg$mean, subseg$p,  pch=19 ,  col="white" , cex= cxcut  - .5 )
      points( subseg$mean, 1-subseg$p,  pch=21 ,  col="blue", lwd=3  , cex=cxcut  )
      points( subseg$mean, 1-subseg$p,  pch=19 , col="white" , cex=cxcut -.5 )
      x <- seq(Sn, 2, .01 ) 
      points( x , ARloh( x , 1 , Sn ) , type='l', col="black"  , lwd=3 )
      points( x , 1- ARloh( x , 1 , Sn ) , type='l', col="black"  , lwd=3 )
    }
    
    return( Sn )
  } else {
    
    if( is.null(cntr)  ){ stop("estimateLOHcurve: missing cntr argument")}
    image( cntr , col=terrain.colors(50) )
    contour( cntr , nlev=50, add=T )
    le<- subseg$end.pos-subseg$start.pos
    cxcut<- as.integer( cut( le, c(100000,1000000,5000000,10000000,20000000,50000000,Inf) ) )/3
    points( x<-subseg$mean, y<-subseg$p,  pch=21 , col="blue", lwd=3 , cex=cxcut  )
    points( subseg$mean, subseg$p,  pch=19 ,  col="white" , cex= cxcut  - .5 )
    points( subseg$mean, 1-subseg$p,  pch=21 ,  col="blue", lwd=3  , cex=cxcut  )
    points( subseg$mean, 1-subseg$p,  pch=19 , col="white" , cex=cxcut -.5 )
    
    cat("Select segments on the LOH curve; right click to finish\n")
    
    selectedPoints<-c()
    zz<-subseg[,c("mean","p")]
    id<-identify( zz[,1],zz[,2] , n=1 , labels="" )
    while( length( id )> 0 ){
      points( zz[id,] , pch=19, col="black") 
      selectedPoints<-rbind( selectedPoints,  data.frame( x=zz[id,1], y=zz[id,2] )  )
      id<-identify(  zz[,1] ,zz[,2] , n=1, labels="" )
    }
    x<-selectedPoints[,1]
    y<-selectedPoints[,2]
    
    nnllss<- nls( y ~ ARloh( x, 1 ,n ) , start=list(n=0.01  ) , upper=list( n=1 ), lower=list(n=0 ) , algo="port" )
    Sn<-  summary(nnllss)$coefficients[1,1]
    x <- seq(Sn, 2, .01 ) 
    points( x , ARloh( x , 1 , Sn ) , type='l', col="black"  , lwd=3 )
    points( x , 1- ARloh( x , 1 , Sn ) , type='l', col="black"  , lwd=3 )
    
    return(Sn)
    
  }
}
    
    