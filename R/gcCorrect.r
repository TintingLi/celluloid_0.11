# my own gc correction


gcCorrect<-function( rangedata , sampletype="normal", span =0.3 , mappability = 0.9, samplesize = 50000 , bprange=6 , maskmap=T ){

  if( sampletype=="tumor"){
    rangedata <- gcCorrect.tumor(  rangedata , span =span , mappability = mappability, samplesize = samplesize , 
                                   bprange=bprange , maskmap=maskmap )
  } else if( sampletype=="normal"){
    rangedata <- gcCorrect.normal(  rangedata , span =span , mappability = mappability, samplesize = samplesize , 
                                    bprange=bprange , maskmap=maskmap )
  } else { stop( "gcCorrect: unknown sampletype (must take one of \"normal\" or \"tumor\" \n") }
  return( rangedata )
}


gcCorrect.tumor<-function( rangedata , span =0.3 , mappability = 0.9, samplesize = 50000 , bprange=6 , maskmap=T ){
      
  # first pass
  cat( "gcCorrect: first pass\n")
  tc<-gcCorrect( rangedata , maskmap=FALSE, sampletype="normal" )
  
  t.seg <- segmentSeqData( tc , k=50 , maskmap=0, skipmeanmap=T )
  
  # will now add a tmp column to tc to include mean segment values 
  
  tmp<-tc[1,]
  binlength<-end( tmp )-start(tmp) + 1
  
  # this extract the segment mean and assigns it to each bin 
  cat(" extracting segment mean and assigning to each bin\n ")
  meaninbins<- apply(   t.seg[,c("chrom","start.pos","end.pos","mean" )] , 1, 
                        function(x){ cat(".")
                                     fr<-floor(as.numeric(x[2])/binlength)*binlength+1
                                     to<-floor(as.numeric(x[3])/binlength)*binlength
                                     if( to>fr ){  # this will exclude segs with length < binlength
                                       tmp<-paste( x[1]  , seq(fr,to,binlength) , sep=":" )
                                       tmp<- data.frame(tmp,as.numeric(x[4]) )
                                       names( tmp)<-c("bin","mean"); return(tmp) }} )
  
  require( data.table) 
  meaninbins <-rbindlist(meaninbins)
  
  # to keep track of the line number in tc
  binlabels <- data.frame( line=1:nrow(tc),
                           bin=paste( tc$space,":", start(tc), sep="" ) )
  
  m<-merge( binlabels , meaninbins, by="bin", all.x=T ) 
  
  # this new entry in tc will hold the segment bins. sm for smooth (from an old code)
  tc$tmp.mean <- NA 
  tc$tmp.mean <- m$mean[ order(m$line) ] 
  
  
  # second pass
  
  cat( "\ngcCorrect: second pass\n")
  tc$gc[ tc$gc<0 ] <- NA 

  # ignore for correction purposes
  
  tc$ignore<-F
  tc$ignore[ is.na(tc$reads.gc) | is.na( tc$gc )] <- T
  tc$ignore[ !is.element( tc$space, paste("chr",1:22,sep="") ) ] <- T 
  
  # ignoring read count that are outliers (in the boxplot sense) for correction purposes  
  
  cat("Identifying and ignoring outlier read counts...")
  
  bp<-boxplot( tc$reads.gc[!tc$ignore] , range=bprange, plot=F )
  tc$ignore[ tc$reads.gc < bp$stats[1,1] | tc$reads.gc > bp$stats[5,1]  ] <- T
  
  tc$ignore[ tc$map < mappability ] <- T 
  
  # random selection of points to speedup loess, done in intervals to help good sampling
  cu<-cut( tc$gc, include.lowest=TRUE, breaks=quantile( tc$gc, seq(0,1,.1),na.rm=T ) )
  keep<-c()
  for( i in levels( cu ) ){
    keep<-c( keep, sample( which( cu==i & !tc$ignore ), floor( samplesize/length(levels(cu)) ) ) )
  }
  
  cat(" done.\nApplying loess correction...")
  
  # Now performing loess by accounting for segment heights
  lo<- loess( tc$reads.gc[keep]/tc$tmp.mean[keep] ~ tc$gc[keep] , span = span )
  
  tmp <- tc$reads.gc/predict(lo, tc$gc)
  tmp[ tmp <0 ]<- NA 
  tmp <- tmp/mean( tmp, na.rm=T )
  tc$reads.gc<-tmp 
  
  tc$copy<-log( tc$reads.gc , 2 )
  
  cat(" done.\nFind the result in $reads.gc\n") 
  
  if( maskmap ){
    tc$ignore <- F 
    tc$ignore[ tc$map < mappability ] <- T 
    tc$reads.gc[ tc$ignore  ]<- NA
    tc$map[ tc$ignore  ]<- NA 
  } else {
    # reverting
    tc$ignore <- F 
  }
  
  return( tc ) 
  
}


gcCorrect.normal<-function( rangedata , span =0.3 , mappability = 0.9, samplesize = 50000 , bprange=6 , maskmap=T ){

   reads <- rangedata$reads

   rangedata$gc[ rangedata$gc<0 ] <- NA 
   reads[ reads<0 ]<-NA 

   # ignore for correction purposes
  
   rangedata$ignore<-F
   rangedata$ignore[ is.na(reads) | is.na( rangedata$gc )] <- T
   rangedata$ignore[ !is.element( rangedata$space, paste("chr",1:22,sep="") ) ] <- T 
   
   # ignoring read count that are outliers (in the boxplot sense) for correction purposes  

   cat("Identifying and ignoring outlier read counts...")

   bp<-boxplot( reads, range=bprange, plot=F )
   rangedata$ignore[ reads < bp$stats[1,1] | reads > bp$stats[5,1]  ] <- T

   rangedata$ignore[ rangedata$map < mappability ] <- T 

   
   # random selection of points to speedup loess, done in intervals to help good sampling
   cu<-cut( rangedata$gc, include.lowest=TRUE, breaks=quantile( rangedata$gc, seq(0,1,.1),na.rm=T ) )
   keep<-c()
   for( i in levels( cu ) ){
     keep<-c( keep, sample( which( cu==i & !rangedata$ignore ), floor( samplesize/length(levels(cu)) ) ) )
   }
   
   cat(" done.\nApplying loess correction...")

   lo<- loess(reads[keep] ~ rangedata$gc[keep], span = span )
  
 if( FALSE ){ 
   plot( rangedata$gc[keep], reads[keep], pch='.' ) 
   or<-order( rangedata$gc[keep] )
   points( rangedata$gc[keep][or], predict(lo,rangedata$gc[keep][or] ), type='l', col="red")  
  }
   rangedata$reads.gc <- rangedata$reads/predict(lo, rangedata$gc)
   rangedata$reads.gc[ rangedata$reads.gc <0 ]<- NA 
   rangedata$reads.gc <- rangedata$reads.gc/mean( rangedata$reads.gc, na.rm=T )

   rangedata$copy<-log( rangedata$reads.gc , 2 )
   cat(" done.\nFind the result in $reads.gc\n") 
 
   if( maskmap ){
     rangedata$ignore <- F 
     rangedata$ignore[ rangedata$map < mappability ] <- T 
     rangedata$reads.gc[ rangedata$ignore  ]<- NA
     rangedata$map[ rangedata$ignore  ]<- NA 
   } else {
     # reverting
     rangedata$ignore <- F 
   }
 
   return( rangedata ) 

}

