# changes: now takes a seg object and a AR object.

# seg has columns chr, start, end, mean 
# chromosomes identifiers are chr1, chr2, .. , chrX, chrY, etc.
# ar dataframe has columns chr, pos, nReadsRef, nReadsAlt
# trange is a range object (IRange) that holds the genomic bins 
# and read counts, etc...

# it has a column called "ignore" which takes T if the bin is 
# to be ignored

prepCopyAr <- function( seg, ar=NULL , tumourrangedata, xonly=FALSE ){

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
 cat("\n merging" )
 m<-merge( binlabels , meaninbins, by="bin", all.x=T ) 
 
 if( FALSE ){ 
   cat( "DEBUG number of lines in m ", nrow( m ) ,"\nDEBUG number of lines in tumourrangedata", dim(tumourrangedata)) 
 }

 cat(".")
 cat("\n adding columns to tc") 
 # this new entry in tumourrangedata will hold the segment bins. sm for smooth (from an old code)
 eval.parent( substitute (
   tumourrangedata$smcopy<- NA ))
 eval.parent( substitute (
   tumourrangedata$smcopy<- m$mean[ order(m$line) ] ))
 eval.parent( substitute (
   tumourrangedata$smcopy[ is.na(tumourrangedata$reads.gc) | tumourrangedata$ignore ] <- NA ))
# need to the same for local to the function?  
 tumourrangedata$smcopy<- NA 
 tumourrangedata$smcopy<- m$mean[ order(m$line) ] 
 tumourrangedata$smcopy[ is.na(tumourrangedata$reads.gc) | tumourrangedata$ignore ] <- NA 

 cat(".")

 # SINGLECELL BRANCH
 if( is.null( ar ) ){
   if( !xonly ){ stop("ar can't be NULL unless xonly=TRUE\n") }
   # creating a fake ar, one het every 10kb 
   cat("\n faking ar")
   ar<-c()
   for( i in 1:nrow( seg ) ) { 
     cat(".")
     pos<- seq( seg[i,4], seg[i,5], 50000 )
     tmp<- data.frame( CHR= seg[i,2],  POS=pos,  REF_COUNT=50,  VAR_COUNT=50 )
     ar<- rbind( ar, tmp ) 
   }
 }

 if( FALSE ) { print( head( tumourrangedata ) )  }
 chpo<- paste( tumourrangedata$space, start(tumourrangedata), sep="-")
 ars <- data.frame( chpo=chpo  ,copy= tumourrangedata$smcopy )
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





 


  
