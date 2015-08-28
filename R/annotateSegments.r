annotateSegments<-function( seg, epp , weights="quadratic" ){
  
  tmpseg<-seg
  tmpseg$labels<-NA
  tmpseg$dist<- NA 
  tmpseg$we<-NA 
  
  #tmp <- apply( epp[, 3:(ncol(epp)-2 )], 1, sum )
  #sel<-tmp==0 | tmp==1 
  
  # distance between "red" lines
  #d<-dist( epp$x[sel] )
  # changed min to mean
  #xdist <-mean( d[d>0] )/2
  # number of red lines expected between two integer values
  sel0 <- apply( epp[, 3:(ncol(epp)-2 )]==0 , 1, all )
  sel1 <- apply( epp[, seq(3,(ncol(epp)-2 ),2)]==1 , 1, all ) & 
                  apply( epp[, seq(4,(ncol(epp)-2 ),2)]==0 , 1, all )
  xdist<- (epp[sel1,"x"]-epp[sel0,"x"])/(nsubcl*2)
  
  tmpseg$xdiff<-xdist
  tmpseg$dist <- NA
  tmpseg$we<-NA 
  
  nseg<-nrow(seg)
  nepp<-nrow(epp)
  
  eoseg<-rbind(as.matrix( seg[,c("mean","p")]  ), as.matrix( epp[,c("x","ar")] ) )
  
  d<-as.matrix(dist(eoseg) )
  d<-d[1:nseg,(nseg+1):(nepp+nseg)] 
  
  mnd<-apply( d, 1, min )
  
  if( weights=="quadratic" ){
    mnd<-weight.quadratic( mnd, xdist=xdist)
  } else if( weights=="linear"){
    mnd<-weight.linear( mnd, xdist=xdist)
  } else {stop("annotateSegments: unknown weights arguments\n")}
  
  
  tmpseg$we<-mnd
   
  # if the distance between a segment and a epp is 0, the function returns 1
  # if it is greated than xdiff (the mid point between two integer value, 
  # it returns 0 
  # I want to calculate some percentage of genome that is captured by epp
  # I use these values as weight
  # linear:
 
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
    m <- epp[sel, ][seq(3, ncol(epp) - 2, 2)]
    p <- epp[sel, ][seq(4, ncol(epp) - 2, 2)]
    lab <- paste(paste(m, p, sep = ""), collapse = "/")
    tmpseg$labels[s]<-lab
  }
  
}

return(tmpseg )

}


