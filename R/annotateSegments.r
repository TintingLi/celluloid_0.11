annotateSegments<-function( seg, epp , weights="quadratic" ){
  
  tmpseg<-seg
  tmpseg$labels<-NA
  tmpseg$dist<- NA 
  tmpseg$we<-NA 
  d<-dist( epp$x )
  xdiff<-min( d[d>0] )/2
  tmpseg$xdiff<-xdiff
  tmpseg$dist <- NA
  tmpseg$we<-NA 
  
  # if the distance between a segment and a epp is 0, the function returns 1
  # if it is greated than xdiff (the mid point between two integer value, 
  # it returns 0 
  # I want to calculate some percentage of genome that is captured by epp
  # I use these values as weight
  # linear:
  if( weights=="linear"){
    we <-function(di){ tmp<- 1-di/xdiff;  return( max( 0, tmp ) ) }
  } else if( weights=="quadratic"){
    we <-function(di){ 
      b<- -2/xdiff 
      a<- -b/(2*xdiff)
      c<- 1
      if( di > xdiff ){ return( 0 )} else {return( a*di^2+b*di+c )}
    }
  } else {stop("annotateSegments: unknown weights arguments\n")}
  
  for( s in 1:nrow( tmpseg ) ){
    
    if( s==1 | s%%10==0 ){
      cat("annotating segment ",s," of ",nrow(tmpseg),"...\n", sep="")
    }
  if( !is.na( tmpseg[s,"mean"] )&!is.na( tmpseg[s,"p"] ) ){
    tmp<- rbind( data.frame( x=tmpseg[s,"mean"], ar=tmpseg[s,"p"] ),
                 epp[,c("x","ar")] )
    d<-as.matrix( dist(tmp) )
    sel<- which( d[1,-1]==min(d[1,-1] ) )[1]
    tmpseg$dist[s]<- d[1,-1][sel][1]
    if( nrow( epp )==6 ){
      tmpseg$we[s]<-we( d[1,-1][sel][1] )
    }
    m <- epp[sel, ][seq(3, ncol(epp) - 2, 2)]
    p <- epp[sel, ][seq(4, ncol(epp) - 2, 2)]
    lab <- paste(paste(m, p, sep = ""), collapse = "/")
    tmpseg$labels[s]<-lab
  }
  
}

return(tmpseg )

}


