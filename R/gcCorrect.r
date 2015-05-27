# my own gc correction

gcCorrect<-function( rangedata , span =0.3 , mappability = 0.9, samplesize = 50000 , bprange=6 ){

   reads <- rangedata$reads

   rangedata$gc[ rangedata$gc<0 ] <- NA 
   reads[ reads<0 ]<-NA 

   # ignore for correction purposes
  
   rangedata$ignore<-F
   rangedata$ignore[ is.na(reads) | is.na( rangedata$gc )] <- T 
   rangedata$ignore[ !is.element( rangedata$space, paste("chr",1:22,sep="") ) ] <- T 
   
   # ignoring read count that are outliers (in the boxplot sense) for correction purposes  

   cat("Identifying and ignoring outlier read counts...")

   bp<-boxplot( reads, rangedata=bprange, plot=F )
   rangedata$ignore[ reads < bp$stats[1,1] | reads > bp$stats[5,1]  ] <- T

   rangedata$ignore[ rangedata$map < mappability ] <- T 

   keep<- which( !rangedata$ignore )
   keep<- sample( keep, samplesize )

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
   return( rangedata ) 

}

