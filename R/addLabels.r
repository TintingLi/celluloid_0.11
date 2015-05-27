
### TODO 
### PROBLEM IF SUBCLONES ARE IN EUQAL FREQ: 

# epp: expected peaks 

addLabels<-function( epp , selectedPoints =NULL ,manual=F, cex=.5, bg="white",...  ){

 selpoints<-c()
 if( !is.null(selectedPoints) ){
  names(selectedPoints)[1:2]<-c("x","ar")

  for( i in 1:nrow(selectedPoints) ){
    # select the point from epp that is closest to selectedPoint
  
    d<-as.matrix( dist( rbind(selectedPoints[i,c("x","ar")], epp[,c("x","ar") ] ) )  )
    sel<-which( d[,1][-1] == min( d[,1][-1] ) )
    m<-epp[sel,][seq(3,ncol(epp)-2,2) ]
    p<-epp[sel,][seq(4,ncol(epp)-2,2) ]
    leg <-  paste( paste(m,p, sep="" ), collapse="/") 
    legend( epp[sel,"x"], epp[sel,"ar"], legend=leg  , cex=cex, bg=bg, ... )
   }
 }
 if( manual ){ 
   cat( "left click on points to add a label; right click to finish\n") 
   sel <-identify( epp$x, epp$ar , n=1 , labels=""  )
   while( length(sel) > 0 ){
     m<-epp[sel,][seq(3,ncol(epp)-2,2) ]
     p<-epp[sel,][seq(4,ncol(epp)-2,2) ]
     leg <-  paste( paste(m,p, sep="" ), collapse="/") 
     legend( epp[sel,"x"], epp[sel,"ar"], legend=leg , xjust=1, cex=cex, bg=bg,... )
     selpoints<-rbind( selpoints, c( epp[sel,"x"], epp[sel,"ar"] ) )
     sel<-identify( epp$x, epp$ar , n=1 , labels="" )
   }
 }
 if( !is.null(nrow(selpoints)) ){
  colnames(selpoints)<-c("x","ar")
  invisible(selpoints)
 }
 
}


