#source("meanMapInSeg.r") 

segmentSeqData<-function( rangedata, gamma=500, kmin=100, maskmap=.8 , maskadj=FALSE  ){

 binlength<- end(rangedata[1,])-start(rangedata[1,])+1 

 tmp<-data.frame(chrom=rangedata$space, pos=start(rangedata), rc=rangedata$reads.gc )

 tmp$chrom<- factor( as.character(tmp$chrom), levels=c(paste("chr",1:22,sep=""),"chrX","chrY" ) )
 tmp$chrom<-as.numeric( tmp$chrom )

 sel<- rangedata$map< maskmap  
 if( maskadj & FALSE ){
  # picking adjacent bins as well; might be a problems when switching chr
  sell<-c( F, sel[-length(sel)] )
  selr<-c( sel[-1], F )
  sel<- sel|sell|selr
 }

 tmp<-tmp[!sel,] 
 tmp<- tmp[ !is.na(tmp$rc),] 

 require(copynumber)

 tmp.win<-winsorize( tmp )
 tmp.seg <- pcf( data=tmp.win, gamma=gamma , kmin=kmin    , digits=4 )
 
 # start and end are based on "probe" positions;
 # I define a probe to be the start of a bin,  extending the end.pos to the end of the bins

 tmp.seg$end.pos<-tmp.seg$end.pos + binlength -1

 tmp.seg$chrom<-paste("chr", tmp.seg$chrom, sep="")   

 tmp.seg.meanmap <-meanMapInSeg(  tmp.seg , rangedata  ) 

 return( tmp.seg.meanmap ) 

}


