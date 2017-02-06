
# copyAr created with prepCopyAr 
#   chr    pos   copy ref alt

showTumourProfile<-function(copyAr, nx=200, ny=50, maxPoints=50000, selected=NULL, 
                            flatten=1, xlim=c(0,2), nlev=50 , xaxt="n", seed=NULL, 
                            nopoints=F , hist2d=F  , chr=NULL , noise=NULL, 
                            col=terrain.colors(50) , onlypoints=F ,plot=T ){

  if( !is.null(seed) ){ set.seed(seed) }
  
  #require(gplots)
  #require(MASS) 
  # SINGLECELL BRANCH
  fakeAR<-NULL ; N<-100
  # if there are no variation in ref and alt, assume that xonly is on
  if( var( copyAr$ref )==0 & var( copyAr$alt )==0 ){
    set.seed(12345) 
    # faking an ar, forcing it to be symmetrical 
    fakeAR<- rep( NA, nrow( copyAr ) )
    # lenght of idx1 will be >= length( idx2 )
    idx1<- seq( 1, nrow( copyAr), 2 )
    idx2<- seq( 2, nrow( copyAr), 2 )
    rb<- rbinom( length( idx1) , N , .5 )
    rb[rb<N/2 ]<-N-rb[rb<N/2 ] 
    fakeAR[ idx1 ]<- rb
    fakeAR[ idx2 ]<- N -rb[ 1:length(idx2) ] 
    copyAr$ref<- fakeAR
    copyAr$alt<- N -copyAr$ref
  }  
  
  
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
      if( is.null( fakeAR ) ) {
        plot(  tmpCopyAr$copy[sa] , tmpCopyAr$ar[sa] , pch='.', 
               col="gray", xlim=xlim , xaxt = xaxt , xlab="Copy number" , ylab="Het AR (%RefAlleles)", ylim=c(0,1)  ) 
      } else {
        plot(  tmpCopyAr$copy[sa] , tmpCopyAr$ar[sa] , pch='.', 
               col="gray", xlim=xlim , xaxt = xaxt , xlab="Copy number" , yaxt="n", ylab="", ylim=c(0,1)  ) 
      }
      
    } else {
      plot.new() 
      plot.window( xlim= xlim , ylim=c(0,1) )
      if( is.null( fakeAR ) ){
        yaxt<-"s"
        title(  xlab="Copy number" , ylab="Het AR (%RefAlleles)"  )
      } else {
        yaxt<-"n" 
        title(  xlab="Copy number"   )
      }
      axis(2) 
      box() 
    }
    
    if( !onlypoints ){
      # cntr$z[ cntr$z < .05* max( cntr$z , na.rm=T) ] <- 0 
      image(  cntr$x, cntr$y, cntr$z ,add=T, col=col  , xaxt="n", yaxt=yaxt )
      contour( cntr$x, cntr$y, cntr$z , xlim=xlim , nlev=nlev  , add=T )
    }
  }
  invisible( cntr )
}





