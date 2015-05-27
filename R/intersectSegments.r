 
intersectSegments<-function( seg1, seg2 ){

  interseg<-c()

  for( ch in union( seg1$chrom , seg2$chrom) ){

   for( arm in c("p","q")){

    sub1<-seg1[ seg1$chrom==ch & seg1$arm==arm,]
    sub2<-seg2[ seg2$chrom==ch & seg2$arm==arm,]

    #names(sub1)[7:ncol(sub1)]<-  paste( names(sub1)[7:ncol(sub1)] , ".1", sep="")
    #names(sub2)[7:ncol(sub2)]<-  paste( names(sub2)[7:ncol(sub2)] , ".2", sep="")

    if( nrow(sub1)+nrow(sub2)>0 ){

     bounds<-sort(unique( c( sub1$start.pos,  sub1$end.pos+1 ,  sub2$start.pos,  sub2$end.pos+1 ) ) )

     for( i in 1:(length(bounds)-1) ){

       s<-bounds[i]
       e<-bounds[i+1]

       sel1<- sub1$start.pos<=s & sub1$end.pos+1>=e 
       sel2<- sub2$start.pos<=s & sub2$end.pos+1>=e 

       if( sum(sel1) ){
        line1<-as.data.frame(sub1[sel1, 7:ncol(sub1) ])
       } else { 
        line1<-as.data.frame( t( rep( NA, ncol(sub1)-6 ) ) )
        }
       if( sum(sel2) ){
        line2<-as.data.frame(sub2[sel2, 7:ncol(sub2) ])
       } else { 
        line2<-as.data.frame( t( rep( NA, ncol(sub2)-6 ) ) )
       }
       colnames(line1)<-names( sub1 )[7:ncol(sub1)] 
       colnames(line2)<-names( sub2 )[7:ncol(sub2)] 

       interseg<-rbind( interseg,  data.frame( sampleID=NA, chrom=ch, arm=arm,
                                        start.pos=s, end.pos=e-1, size=e-s,
                                        line1, line2 ) )

    }
   }
  }
 }
return(interseg)
}


