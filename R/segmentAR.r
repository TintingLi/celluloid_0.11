# ar is data.frame with columns 
#   CHR   POS REF_COUNT VAR_COUNT
  
# bin penalty by default; this segmenting is to identify segments that alternate between LOH and normal 
segmentAR<-function( ar,  tumourrangedata=NULL, gamma = 500, kmin=100,  maskmap=.8  ){

  ar.tmp<-ar

  # excluding hets that belong in a masked bin
  if( maskmap>0 & !is.null( tumourrangedata ) ){
    binlength<- end( tumourrangedata[1,] ) -start( tumourrangedata[1,] ) +1 
    # in which bin does the het belong to 
    bin<-paste( ar[,1], binlength*floor( ar[,2]/binlength )+1, sep="-") 
    sel<- tumourrangedata$map < maskmap 
    mask <- paste( tumourrangedata$space, start(tumourrangedata), sep="-")[sel]
    ar.tmp<-ar[!is.element( bin, mask ),] 
  }
 

  tmp<-data.frame(chrom=ar.tmp[,1], pos=ar.tmp[,2] , ar =ar.tmp[,3]/(ar.tmp[,3]+ar.tmp[,4]) ) 
  # selecting autosomal chromosomes
  sel<-is.element( tmp$chrom, paste("chr",1:22,sep="") )
  tmp<-tmp[sel,]
  # relabel chromosomes into integers
  tmp$chrom<- factor( as.character(tmp$chrom), levels=c( paste("chr",1:22,sep="") ) )
  tmp$chrom<-as.numeric( tmp$chrom )
  # don't know why this was there...
  tmp<- tmp[ !is.na( tmp$chrom ),] 
  tmp<-tmp[ order( tmp$chrom, tmp$pos),] 
  
  # recomputing AR so that all are <= 0.5 
  sel<-tmp$ar>.5
  tmp$ar[sel]<-1-tmp$ar[sel]
  
  #require(copynumber)
  tmp.win<-winsorize( tmp )
  ar.seg<- pcf( data=tmp.win, gamma=5000 , kmin=100    , digits=4 )
  ar.seg$chrom<-paste("chr", ar.seg$chrom, sep="")   

  colnames(ar.seg)[7]<-"meanar"
  
 return( ar.seg )


}

