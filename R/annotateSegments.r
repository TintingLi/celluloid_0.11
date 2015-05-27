annotateSegments<-function( seg, epp ){

 tmpseg<-seg
 tmpseg$labels<-NA

 for( s in 1:nrow( tmpseg ) ){

 if( s==1 | s%%10==0 ){
  cat("annotating segment ",s," of ",nrow(tmpseg),"...\n", sep="")
 }
 if( !is.na( tmpseg[s,"mean"] )&!is.na( tmpseg[s,"p"] ) ){
  tmp<- rbind( data.frame( x=tmpseg[s,"mean"], ar=tmpseg[s,"p"] ),
          epp[,c("x","ar")] )
  d<-as.matrix( dist(tmp) )
  sel<- which( d[1,-1]==min(d[1,-1] ) )[1]
   m <- epp[sel, ][seq(3, ncol(epp) - 2, 2)]
   p <- epp[sel, ][seq(4, ncol(epp) - 2, 2)]
   lab <- paste(paste(m, p, sep = ""), collapse = "/")
   tmpseg$labels[s]<-lab
 }

 }

 return(tmpseg )

}


