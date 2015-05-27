# changes: now takes a seg object and a AR object.

# seg has columns chr, start, end, mean 
# chromosomes identifiers are chr1, chr2, .. , chrX, chrY, etc.
# ar dataframe has columns chr, pos, nReadsRef, nReadsAlt
# trange is a range object (IRange) that holds the genomic bins 
# and read counts, etc...

# it has a column called "ignore" which takes T if the bin is 
# to be ignored

prepCopyAr <- function( seg, ar, trange ){
 require( data.table) 

 tmp<-trange[1,]
 binlength<-end( tmp )-start(tmp) + 1

 # this extract the segment mean and assigns it to each bin 
 meaninbins<- apply( data.frame( seg ), 1, 
                      function(x){ 
                        fr<-floor(as.numeric(x[2])/binlength)*binlength+1
                        to<-floor(as.numeric(x[3])/binlength)*binlength
                        tmp<-paste( x[1]  , seq(fr,to,binlength) , sep=":" )
                        tmp<- data.frame(tmp,as.numeric(x[4]) )
                        names( tmp)<-c("bin","mean"); return(tmp) } )
 meaninbins <-rbindlist(meaninbins)

 # to keep track of the line number in trange

 binlabels <- data.frame( line=1:nrow(trange), bin=paste( trange$space,":", start(trange), sep="" ) )
 m<-merge( binlabels , meaninbins, by="bin", all.x=T ) 
 
 # this new entry in trange will hold the segment bins. sm for smooth (from an old code)
 trange$smcopy<- NA 
 trange$smcopy<- m$mean[ order(m$line) ] 

 trange$smcopy[ trange$ignore ] <- NA 

 ars <- data.frame( chpo=paste( trange$space, start(trange), sep="-") ,copy= trange$smcopy )
 ars<-ars[ !tc$ignore,] 

 alratio<-ar
 alratio<-alratio[ alratio[,3]+alratio[,4]>0 ,] 

 alratio$chpo <- paste( alratio[,1] , floor( alratio[,2]/binlength )*binlength +1  , sep="-") 

 copyAr <-merge( copyAr, alratio , by="chpo" )

 copyAr<-copyAr[ , c(3,4,2,5,6)] 
 names(copyAr)<-c("chr","pos","copy","ref","alt") 
 copyAr<-copyAr[order(copyAr$chr, copyAr$pos),] 
 copyAr<-copyAr[ !is.na(copyAr$copy),] 

 rownames(copyAr)<-1:nrow(copyAr)
 
 return(copyAr) 

}



 


  