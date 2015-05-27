
showSTumourProfile<-function(scopyAr, nx=200, ny=50, maxPoints=50000, selected=NULL, 
                            flatten=.25, xlim=c(0,2), nlev=50  ){
 #require(MASS)
 if( is.finite(maxPoints) ){
   sa<-sample( 1:nrow(scopyAr), min( maxPoints, nrow(scopyAr)) , replace=F )
 } else { 
   sa<-1:nrow( scopyAr )
 }
 
 cntr <-kde2d( scopyAr[sa,"copy"] , scopyAr[sa,"sar" ]  , n=c( nx, ny ) , lims=c( xlim, 0,1 ) )
 plot(  scopyAr[sa,1] , scopyAr[sa,2] , pch='.', col="gray", ylim=c(0,1), 
         xlim=xlim , xaxt = "n" , xlab="Copy number" , ylab="SV AR (%RefAlleles)" ) 
 image(  cntr$x, cntr$y, cntr$z^flatten ,add=T , col=terrain.colors(50) )
 contour( cntr$x, cntr$y, cntr$z^flatten , xlim=xlim , nlev=nlev  , add=T )

}
