
# copyAr created with prepCopyAr 
#   chr    pos   copy ref alt

showTumourProfile<-function(copyAr, nx=200, ny=50, maxPoints=50000, selected=NULL, 
                            flatten=1, xlim=c(0,2), nlev=50 , xaxt="n", seed=NULL, 
                            nopoints=F , hist2d=F  , chr=NULL , noise=NULL, 
                            col=terrain.colors(50) , onlypoints=F ,plot=T ){
  
  #require(gplots)
  #require(MASS) 
  
  if( is.null(chr) ){
    tmpCopyAr<-copyAr
  } else {
    sel<-copyAr$chr==chr
    tmpCopyAr<-copyAr[sel,]
  }
  tmpCopyAr$ar<-tmpCopyAr$ref/( tmpCopyAr$ref+tmpCopyAr$alt )
  
  tmpCopyAr<-tmpCopyAr[ !is.na( tmpCopyAr$copy ) , ]
  
  if( !is.null( noise ) ){
    tmpCopyAr$copy <- tmpCopyAr$copy + rnorm( length(tmpCopyAr$copy), 0, noise) 
  }
  
  if( !is.null(seed) ){ set.seed(seed) }
  
  sel<-tmpCopyAr$copy >= xlim[1] & tmpCopyAr$copy <= xlim[2]
  
  #require(MASS)
  if( is.finite(maxPoints) ){
    sa<-sample( (1:nrow(tmpCopyAr))[sel]  , min( maxPoints, nrow( tmpCopyAr[sel,]   )) , replace=F )
  } else { 
    sa<-(1:nrow( tmpCopyAr ))[sel] 
  }
  
  if( !hist2d ){
    cntr <-kde2d( tmpCopyAr$copy[sa] , tmpCopyAr$ar[sa]  , n=c( nx, ny ) , lims=c( xlim, 0,1 ) )
  } else {
    cntr <-hist2d( tmpCopyAr$copy[sa] , tmpCopyAr$ar[sa]  , nbins=c( nx, ny ) , show=F )
    cntr$z<-cntr$counts 
  }
  
  cntr$z<-cntr$z^flatten 
  
  if( onlypoints ){ nopoints<- F }
  if( plot ){
    if( !nopoints ){
      # add points; for use with identify() 
      plot(  tmpCopyAr$copy[sa] , tmpCopyAr$ar[sa] , pch='.', 
             col="gray", xlim=xlim , xaxt = xaxt , xlab="Copy number" , ylab="Het AR (%RefAlleles)", ylim=c(0,1)  ) 
    } else {
      plot.new() 
      plot.window( xlim= xlim , ylim=c(0,1) )
      title(  xlab="Copy number" , ylab="Het AR (%RefAlleles)"  )
      axis(2) 
      box() 
    }
    
    
    if( !onlypoints ){
      # cntr$z[ cntr$z < .05* max( cntr$z , na.rm=T) ] <- 0 
      image(  cntr$x, cntr$y, cntr$z ,add=T, col=col  , xaxt="n" )
      contour( cntr$x, cntr$y, cntr$z , xlim=xlim , nlev=nlev  , add=T )
    }
  }
  invisible( cntr )
}
