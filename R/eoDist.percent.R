
peakProximity.percent <-function( v, segments , verbose=T, Sn=NULL , ... ){
  
  if( is.null(Sn) ){
    S<-v[1] 
    t<- tail(v, -1 ) 
    t<- c( t, 1-sum(t) )
  } else {
    t<-v
    if( t[1]<0.001 ){ t[1]<-0.001 }
    t<-c( t, 1-sum(t) )
    S<-Sn/t[1]
  }
  
  totDist <- eoDist.percent( segments , S=S, t=t , ... ) 
  
  if( totDist <  thisMn ){
    thisMn <<- totDist 
  }
  
  # IN ALL OUTPUT I USE PERCENT CAPTURE, BUT THE CODE 
  # MINIMIZED 1-PERCENTCAPTURE
  # output 
  ta<- c( 1-totDist, S , t )  
  if( totDist < 1  ){ 
    tmpglobal<-paramSpace 
    paramSpace <<-rbind( tmpglobal, ta ) 
  }
  
  if( DEBUG ){
    cat(verbose,"\n")
    cat(totDist,"\n")
    cat(thisMn,"\n")
  }
  
  if(  (verbose &  totDist<=thisMn )  ){     
    cat(paste( paste("value=",1-totDist, sep="" ) , paste("S=",S, sep=""), 
               paste( "t=c(", paste( t, collapse=", "), ")" ), sep="; "),"\n")
  } 
  return( totDist )
  
}


##################################################################
# here segments is an obsect such as t.ar.seg. It must have columns named mean, p and size
# this only works for one clone models 

eoDist.percent<-function( segments ,  S, t , quadratic=TRUE ,...  ){
  
  if( sum( is.element( c("mean","p","size"), colnames(segments) ) )!=3){
    stop("eoDist.percent: segments must have columns named mean, p and size")
  }
  
  seg<-segments[ !apply( is.na(segments[,c("mean","p") ] ), 1, any ), ]
  
  # first element of t is the % of normal cells 
  nsubcl <- length(t)-1 
  
  if( nsubcl > 1 & FALSE  ){ stop("eoDist.percent only implemented for one clone models") }
  
  if( abs( sum(t) -1 )>1e-8 ){ stop("parameter vector does not sum to 1 (eoDist.percent)") }
  
  # in case % normal cells is 0 
  if( t[1]==0 ){ t[1]<-1e-4; t<-t/sum(t) } 
  if( t[1]==1 ){ t[1]<-.9999; t<-t/sum(t) }
  # I am only allowing decreasing frequencies ( t[2]>=t[3]>=... ) 
  # all % must be positive
  if( nsubcl>1 ){
    if(  all( t>=0 ) & all(t<=1)  & ( any( t[2:(nsubcl)]< t[3:(nsubcl+1)] )   )  ){ return( 1+max( t[3:(nsubcl+1)] ) ) }
  } 
  
  if( any( t<=0 ) ){ return( 1 - min(t) ) }
  if( any( t>=1 ) ){ return( max( t) ) }
  if( S> max( seg[,"mean"] , na.rm=T ) ){ return( 1+S ) }  
  
  epp<-ePeakPos( S=S, t=t, cn=cn, ... )
  
  #tmp <- apply( epp[, 3:(ncol(epp)-2 )], 1, sum )
  #sel<-tmp==0 | tmp==1 
  
  # distance between "red" lines
  #d<-dist( epp$x[sel] )
  # changed min to mean
  #xdist <-mean( d[d>0] )/2
  # xdist has changed, see below
  
  nseg<-nrow(seg)
  nepp<-nrow(epp)
  
  eoseg<-rbind(as.matrix( seg[,c("mean","p")]  ), as.matrix( epp[,c("x","ar")] ) )
  
 
  if( FALSE ){
  # for each segment, I want to determine the closest "red" lines 
  # on both sides, xdist will be calculated in between them
  d<-outer( as.vector(eoseg[,1]), as.vector(eoseg[,1]) , '-')
  d<-d[1:nseg,(nseg+1):(nepp+nseg)] 
  whright<- apply( d, 1, function(v){
                          if( !any( v<0 ) ){
                             return(NA)
                          } else { 
                            return(which( v==max( v[v<0] ) )[1] ) } 
                          } )       
  whleft<- apply( d, 1, function(v){
                          if( !any( v>=0 ) ){
                            return(NA)
                          } else { 
                            return(which( v==min( v[v>=0] ) )[1] ) } 
                          } )
  
  
  xdist<-(epp[whright,"x"]-epp[whleft,"x"] )/2
  xdist[ is.na( xdist ) ] <- min( abs( d[d!=0] ), na.rm=T )
  } 
  if( TRUE ){
    sel0 <- apply( epp[, 3:(ncol(epp)-2 )]==0 , 1, all )
    if( ncol(epp)-2 > 4 ){ # if there are subclones
      sel1 <- apply( epp[, seq(3,(ncol(epp)-2 ),2)]==1 , 1, all ) & 
        apply( epp[, seq(4,(ncol(epp)-2 ),2)]==0 , 1, all )
    } else {
      sel1<- epp[,3]==1 & epp[,4]==0 
    }
    # the x distance between 0 copy and 1 copy
    xdist<- (epp[sel1,"x"]-epp[sel0,"x"])/(nsubcl*2)
  }
  
  d<-as.matrix(stats::dist(eoseg) )
  # distances between segments and epp (only on the x-axis)
  d<-d[1:nseg,(nseg+1):(nepp+nseg)] 
  
  # for each segment, take the min one
  mnd<-apply( d, 1, min )
  
  if( quadratic ){
    mnd<-weight.quadratic( mnd, xdist=xdist)
  } else {
    mnd<-weight.linear( mnd, xdist=xdist)
  }
  if( nsubcl==1 ){ 
  glength<-sum( seg$size )
  percent<-sum( mnd * seg$size )/glength
  } else {
   # here I am trying to capture individual segments, instead of their sizes. 
    # this is because large noisy segments around a peak can pull too much,
    # as a result the second subclone will be in small frequency to capture the noise.
    glength<-nrow(seg)
    percent<-sum( mnd  )/glength
    
    
  }
    
  # so that the function can be minimized 
  return( 1 - percent )
  
}


################################


weight.linear <-function(di, xdist){ 
  we <- 1-di/xdist;
  we[ di>xdist ]<-0
  return(we)
}

weight.quadratic <-function(di, xdist ){ 
  b<- -2/xdist 
  a<- -b/(2*xdist)
  c<- 1
  we<- a*di^2+b*di+c  
  we[ di > xdist ] <- 0 
  return(we)
}
