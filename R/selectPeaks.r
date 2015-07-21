
filterSym<-function(sp, filtersymm = T){
  
  if( !filtersymm  | nrow( sp )== 1 ){ return(sp) }
  tmpsp<-sp
  # first remove points that are close one another
  d <-as.matrix( dist( tmpsp) ) 
  rem<-c()
  for( i in 1:( nrow( tmpsp )-1 )   ){
    # removing the mirror of the point, so that single peaks at .5 are not removed 
    sel<-(i+1):nrow(tmpsp)
    # checking for close solutions
    if( min( d[i,sel ] )  < 0.01  ){
      # i am removing i as opposed to wh, so that I can ignore i in next for
      rem<-c(rem,i)
      wh<- which( d[i, 1:nrow(tmpsp) ] == min(  d[i,sel ] ) ) 
      # taking the average of point i and point closest to i. 
      tmpsp[ wh , 1 ]<- (tmpsp[ wh  ,1 ]+tmpsp[ i  ,1 ]) /2
      tmpsp[ wh  ,2 ]<- (tmpsp[ wh  ,2 ]+tmpsp[ i  ,2 ]) /2
    }
  }  
    tmpsp<- tmpsp[ !is.element( 1:nrow( tmpsp ), rem ), ] 
  if( nrow(tmpsp)==1 ){ 
    return(tmpsp) 
  }
  # next remove points that are symmetrical. 
  tmp<-tmpsp
  tmp$y<-1-tmp$y
  tmpsp<-rbind( tmpsp, tmp )
  rem<-c()
  d <-as.matrix( dist( tmpsp ) ) 
  for( i in 1:( nrow( tmpsp )/2  -1 )   ){
    sel<- ( nrow( tmpsp )/2+i+1):nrow(tmpsp) 
    if( min( d[i,sel ] )  < 0.025  ){
      rem<-c(rem,i)
      wh<- sel[ which( d[i, sel ] == min(  d[i,sel ] ) )  ]
      # taking the average of point i and point closest to i. 
      tmpsp[ wh  ,1 ]<- (tmpsp[ wh  ,1 ]+tmpsp[ i  ,1 ]) /2
      tmpsp[ wh  ,2 ]<- (tmpsp[ wh  ,2 ]+tmpsp[ i  ,2 ]) /2
    }  
  }
  
  if( length( rem ) > 0 ){
    cat( "Filtered out these mirrored points:\n") 
    print( tmpsp[rem,] )
    # points( tmpsp[ rem, 1:2], pch="X", cex=1.5, col="black" )
  }
  rem<-c(rem, ( nrow( tmpsp )/2+1):nrow(tmpsp)  )
  tmpsp<- tmpsp[ !is.element( 1:nrow( tmpsp ), rem ), ]
  return( tmpsp )
}


selectPeaks<-function( cntr  , copyAr , manual=T, getLocalMax=F, percentMax=.05 , nrand=100 , filtersymm=T,
                        autocol="red", autopch=19, autocex=1, manucol="black", manupch=19, manucex=1 ){
  
  #require(fields) 
  cntr.interpol<- function( xy ){ -1*interp.surface( cntr , matrix( xy , nrow=1 ) )  }
  rax<-range(cntr$x)
  ray<-range(cntr$y)
  
  selectedPoints<-c()
  
  if( getLocalMax ){
    q<-percentMax*max( cntr$z )
    x<-runif( 1, rax[1], rax[2] )
    y<-runif( 1, ray[1], ray[2] )
    op<- optim( c(x,y) , cntr.interpol )
    localMax <- c( op$par, op$value ) 
    for( i in 2:nrand ){ 
      #x<-runif( 1, rax[1], rax[2] )
      #y<-runif( 1, ray[1], ray[2] )
      x<-sample( cntr$x, 1, prob=apply( cntr$z, 1, sum ) )
      y<-sample( cntr$y, 1, prob=apply( cntr$z, 2, sum ) )
      op<- optim( c(x,y) , cntr.interpol )
      lmx <- c( op$par, op$value ) 
      di <- as.matrix( dist( rbind( lmx, localMax ) ) )[,1][-1]
      if( min(di) > 1e-6 ){
        localMax<-rbind( localMax, lmx )
      }
    }
    localMax<-localMax[ localMax[,3] < -q , 1:2 ] 
    selectedPoints<-rbind( selectedPoints,data.frame( x=localMax[,1], y=localMax[,2])  )
    selectedPoints<-filterSym( selectedPoints , filtersymm  )
    points( selectedPoints,  col=autocol, pch=autopch, cex=autocex ) 
    
  }
  
  if( manual ){
    # added this to force a finer grid and not be dependent on points in copyAr
    rx<-range( copyAr[,1] )
    yy<-data.frame(by=1,seq( 0,1, len=101 ) )
    xx<-data.frame(by=1,seq( rx[1], rx[2], len=501 ) )
    zz<-merge( xx,yy, by="by")[,2:3]
    names(zz)<-c("copy","ar")
    cat( "left click on peaks, right click to finish\n") 
    id<-identify( zz[,1] ,zz[,2] , n=1 , labels="" )
    while( length( id )> 0 ){
      points( zz[id,] , pch=19, col="black") 
      selectedPoints<-rbind( selectedPoints,  data.frame( x=zz[id,1], y=zz[id,2] )  )
      id<-identify(  zz[,1] ,zz[,2] , n=1, labels="" )
    }
  }
  
  # removing points that are symmetrical 
  ##      selectedPoints <- sp   
   selectedPoints<-filterSym( selectedPoints, filtersymm  ) 
  
  if( !is.null( nrow(selectedPoints) ) ){
    selectedPoints<-selectedPoints[ order( selectedPoints[,1] ),]
  }
  return(selectedPoints)
  
}

##########################################################################

removePeaks<-function(selectedPoints){
  
  cat( "left click on peaks, right click to finish\n") 
  id<-identify( selectedPoints[,1] ,selectedPoints[,2] , n=1 , labels="" )
  while( length( id )> 0 ){
    points( selectedPoints[id,] , pch="X", col="black") 
    remove<-c(remove, id)
    id<-identify(  selectedPoints[,1] ,selectedPoints[,2] , n=1, labels="" )
  }
  return( selectedPoints[ !is.element( 1:nrow(selectedPoints), remove) , ] )
}

##########################################################################



