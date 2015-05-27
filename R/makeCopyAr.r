
makeCopyAr<-function( seg, ar, tumourrangedata=NULL , maskmap=.8 ){

 # RC in bin versus AR

  binlength<-NA

  ar.tmp<-ar 
  ar.tmp$ar<-ar[,3]/(ar[,3]+ar[,4]) 

  if( maskmap>0 & !is.null( tumourrangedata) ){
    binlength<- end( tumourrangedata[1,] ) -start( tumourrangedata[1,] ) +1 
    # in which bin does the het belong to 
    bin<-paste( ar.tmp[,1], binlength*floor( ar.tmp[,2]/binlength )+1, sep="-") 
    sel<- tumourrangedata$map < maskmap 
    mask <- paste(tumourrangedata$space, start(tumourrangedata), sep="-")[sel]
    ar.tmp<-data.frame( bin=bin, ar.tmp )
    ar.tmp<-ar.tmp[!is.element( bin, mask ),] 
  } else {
    binlength<-1000
    bin<-paste( ar.tmp[,1], binlength*floor( ar.tmp[,2]/binlength )+1, sep="-") 
    ar.tmp<-data.frame( bin=bin, ar.tmp  )
  }

    # this extract the segment mean and assigns it to each bin 
  meaninbins<- apply( data.frame( seg[,c(2,4,5,7)] ), 1, 
                      function(x){ 
                          tmp<-paste( x[1]  , seq( as.numeric(x[2]), as.numeric(x[3])-1, binlength) , sep="-" ); 
                          tmp<- data.frame(tmp,as.numeric(x[4]) ); names( tmp)<-c("bin","mean"); return(tmp) } )

  
  require( data.table) 
  meaninbins <-rbindlist(meaninbins)
  setnames(meaninbins,c("bin", "copy"))
  

  rcar<- merge( ar.tmp, meaninbins, by="bin")
  rcar <- rcar[ rcar[,5]+rcar[,6]>0 ,] 

  copyAr<-list( copyAr=rcar[, c("copy","ar","CHR","POS")] )
  
  return( copyAr )

}


