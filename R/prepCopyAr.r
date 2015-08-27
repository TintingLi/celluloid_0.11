# changes: now takes a seg object and a AR object.

# seg has columns chr, start, end, mean 
# chromosomes identifiers are chr1, chr2, .. , chrX, chrY, etc.
# ar dataframe has columns chr, pos, nReadsRef, nReadsAlt
# trange is a range object (IRange) that holds the genomic bins 
# and read counts, etc...

# it has a column called "ignore" which takes T if the bin is 
# to be ignored

prepCopyAr <- function( seg, ar, tumourrangedata ){

 #require( data.table) 

 tmp<-tumourrangedata[1,]
 binlength<-end( tmp )-start(tmp) + 1

 # this extract the segment mean and assigns it to each bin 
 cat(" extracting segment mean and assigning to each bin\n ")
 meaninbins<- apply(   seg[,c("chrom","start.pos","end.pos","mean" )] , 1, 
                      function(x){ cat(".")
                        fr<-floor(as.numeric(x[2])/binlength)*binlength+1
                        to<-floor(as.numeric(x[3])/binlength)*binlength
                        if( to>fr ){  # this will exclude segs with length < binlength
                         tmp<-paste( x[1]  , seq(fr,to,binlength) , sep=":" )
                         tmp<- data.frame(tmp,as.numeric(x[4]) )
                         names( tmp)<-c("bin","mean"); return(tmp) }} )
 cat(".")
 meaninbins <-rbindlist(meaninbins)

 # to keep track of the line number in tumourrangedata
 cat(".")
 binlabels <- data.frame( line=1:nrow(tumourrangedata),
                          bin=paste( tumourrangedata$space,":", start(tumourrangedata), sep="" ) )
 cat(".")
 m<-merge( binlabels , meaninbins, by="bin", all.x=T ) 
 
 # this new entry in tumourrangedata will hold the segment bins. sm for smooth (from an old code)
 tumourrangedata$smcopy<- NA 
 cat(".")
 tumourrangedata$smcopy<- m$mean[ order(m$line) ] 

 cat(".")
 tumourrangedata$smcopy[ is.na(tumourrangedata$reads.gc) | tumourrangedata$ignore ] <- NA 
 cat(".")

 ars <- data.frame( chpo=paste( tumourrangedata$space, start(tumourrangedata), sep="-") ,copy= tumourrangedata$smcopy )
 ars<-ars[ !is.na(tumourrangedata$reads.gc) & !tumourrangedata$ignore,] 
 cat(".")

 alratio<-ar
 alratio<-alratio[ alratio[,3]+alratio[,4]>0 ,] 
 cat(".")

 alratio$chpo <- paste( alratio[,1] , floor( alratio[,2]/binlength )*binlength +1  , sep="-") 

 cat("\n merging with ref/alt read counts")
 copyAr <-merge( ars, alratio , by="chpo" )
 copyAr$ar<- copyAr[,5]/(copyAr[,5]+copyAr[,6])
 
 copyAr<-copyAr[ , c(2,7,3,4,5,6)] 
 names(copyAr)<-c("copy","ar","chr","pos","ref","alt") 
 copyAr<-copyAr[order(copyAr$chr, copyAr$pos),] 
 copyAr<-copyAr[ !is.na(copyAr$copy),] 

 rownames(copyAr)<-1:nrow(copyAr)
 cat(" ...done.\n")
 return(copyAr) 

}





 


  
