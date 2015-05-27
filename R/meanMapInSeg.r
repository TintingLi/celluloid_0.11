
meanMapInSeg<-function( seg, rangedata ){
 
 seg$meanmap<-NA
 s<-start(rangedata)
 e<-end(rangedata)
 c<-rangedata$space 
 m<-rangedata$map 
 for( i in 1:nrow( seg) ){
  if( i%%10==0 | i==1 )
    cat("calculating mean mappability in segment",i,"of",nrow(seg),"\n")
  sel<- c == seg[i,"chrom"] & s>=seg[i,"start.pos"] & e<= seg[i,"end.pos"] 
  seg[i,"meanmap"]<-mean( m[sel], na.rm=T )
 }
 return(seg)
}

