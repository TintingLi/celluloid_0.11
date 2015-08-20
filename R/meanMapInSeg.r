
meanMapInSeg<-function( seg, rangedata ){
 
 seg$meanmap<-NA
 s<-start(rangedata)
 e<-end(rangedata)
 c<-rangedata$space 
 m<-rangedata$map 
 for( i in 1:nrow( seg) ){
  if( i%%25==0 | i==1 )
    cat("calculating mean mappability in segment",i,"of",nrow(seg),"\n")
  chr<-seg[i,"chrom"]
  if( chr=="chr23" ){ chr="chrX" } 
  if( chr=="chr24" ){ chr="chrY" }
  sel<- c == chr & s>=seg[i,"start.pos"] & e<= seg[i,"end.pos"] 
  seg[i,"meanmap"]<-mean( m[sel], na.rm=T )
 }
 return(seg)
}

