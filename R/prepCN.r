#########################################################
# This function creates a matrix (cn) of all copy number possibilities in normal cells and subclones.
# It creates a global dataframe called cn.
 
# only those configurations are retained: lines with a maximum copy 
# number of maxc (default 2 ) and
# such that there are no more than maxsubcldiff difference in copy 
# number between any two subclones.

prepCN<-function( maxc , nsubcl=1 , maxsubcldiff=NULL  , maxlines=10000  ){

 if( maxc^nsubcl > maxlines ){ stop( "maxc^nsubcl too large") }

 cn<- expand.grid( rep( list( 0:maxc ), nsubcl ) )
 names(cn)<-paste( "T", 1:nsubcl , sep="")

 if( nsubcl>1 ){
  if( maxsubcldiff>=1 ){
   sel<- apply( cn,1, function(x){ max( dist( x ) ) } ) <= maxsubcldiff 
   cn<-cn[sel,] 
  } else { 
    # if maxsubcl<1 
    # one subclone can have a ratio of copy number of no more than 1/maxsubcldiff
    # relative the copies of the other; 
    # ( when maxsubcldiff < 1)
    # exception made for 0 and 1 copies in one, that are ignored 
    sel<-apply( cn,1, function(x){ mn<-max(1,min(x)) ; mx<-max(x); return( mn/mx>=maxsubcldiff ) } ) 
    cn<-cn[sel,]
  }
 }

  # In normal (N), all segments are assumed to have two copies 
  cn <- data.frame( N=2, cn )
  cnw <<-rep(1, nrow(cn))

  cn<<- cn

  invisible(cn)

}




