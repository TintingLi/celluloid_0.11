makeTransparent<-function(someColor, alpha=0)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

##  tlwd=5; tlty=1;tcol=NULL; nlwd=3; nlty=3
plotSegment <- function( tumourrangedata, segments, ar=NULL, n.rc.seg =NULL , columns=NULL, maskmap=NULL, 
                       file="Rplot%03d.pdf", title=NULL, chr=NULL, perpage=4 , 
                       layoutmat=NULL, width= 8.5, height=11 , 
                       ylim=c(-.5,8 ), normal=F, 
                       tlwd=5, tlty=1, tcol=NULL, nlwd=3, nlty=3, ncol=gray(.5) , annotation=NULL, cex.annotation=1 ,...){


if( !is.null(file) ){
   pdf( file=file, width=width, height=height, ... )
}


#   column<-colnames(tumourrangedata)==column

if( is.null(columns) ){
  columns=c("icopy", "imean" )
}
 if( length( columns==2 ) ){
   tccolumn<-colnames(tumourrangedata)==columns[1] 
   segcolumn<-colnames(segments)==columns[2]
 } else {stop("columns argument must be of length 2\n")}

if( sum(tccolumn)!=1 ){ stop("could not find column named", columns[1], "in argument tumourrangedata\n") }
if( sum(segcolumn)!=1 ){ stop("could not find column named", columns[2], "in argument segments\n") }


   chrlist<- c( paste( "chr",1:22, sep="") , "chrX", "chrY" ) 
   if( !is.null( chr ) ){
    chrlist <- chr 
   }

 if( is.null(layoutmat) ){
   tmp<-c()
   for( i in 1:perpage ){
     tmp<-c( tmp, c( 2*i-1, 2*i-1, 2*i-1, 2*i ) )
   }
   layout( matrix( tmp, ncol=1 ) )
 } else {
   layout( layoutmat )
 }

   if( is.null( tcol ) ){ 
#     ra<-c( "red","red", "black", rep("blue", 7 ) )
     rb<-rainbow(24)
     ra<-c( rb[1], rb[4], "black", rb[6], rb[9], rb[13], rb[15], rb[17], rb[19] ,rb[20] )
   } else if ( length( tcol ) == 12 ){
     ra<- tcol 
   } else if ( length( tcol ) == 1 ){
     ra<-rep(tcol, 12 )
   } else {
     stop("tcol has to be NULL or have length of 1 or 12")
   }

   segcol<-ra[as.integer(cut( segments[,segcolumn] , br=c( -Inf, seq( 0.5,8.5,1 ), Inf ) ) )] 
   chrlist<-intersect( chrlist, unique( segments$chrom  ) ) 
   for( chr in chrlist ){

    cat( "plotting" , chr,"\n")

       sel<- tumourrangedata$space==chr 
       x<-( start(  tumourrangedata[sel,] )+ end(  tumourrangedata[sel,] ) )/2
       XLIM<-c(0,max(x))
       y<- tumourrangedata[sel, tccolumn][[1]] 

       if( !is.null( maskmap ) ){
         mask<-tumourrangedata$map[sel]<maskmap
         x<-x[!mask]
         y<-y[!mask]
       }
    
       par( mar=c(0.5,5,6,1 ) )
       plot.new()
       plot.window( ylim=ylim , xlim=XLIM ) 
       title( ylab="copy number" , xlab="" )
       box()
       axis(2)
       # plot( x, y, pch='.', col="gray" , ylim=ylim , xlim=XLIM, xlab="", xaxt="n" , ylab="copy number"  , xaxs="i" )

       main<-chr
       if( !is.null( title ) ){ main<-paste( title,chr, sep="/" ) }
       title( main=main  , line=4.5)
       abline( h=0:ceiling(max(ylim)), lty=1 , col="black" )
       sel<- segments$chrom ==chr 
       if( sum(sel)>0 ){
        subseg<-segments[sel,]

        for( s in 1:nrow(subseg)){
            selx<-x>= subseg[s,"start.pos"] & x<= subseg[s,"end.pos"]
            points( x[selx], y[selx], pch='.', col=makeTransparent(segcol[sel][s], alpha=50) )
            shadow<-"black"
            if( segcol[sel][s]=="black" ){ shadow<-"white" }
            if( !is.null(maskmap) )
              if( subseg[s,"meanmap"]>maskmap ){
                lines( c( subseg[s,"start.pos"], subseg[s,"end.pos"] ), c( subseg[s,segcolumn], subseg[s,segcolumn] ) , 
                   lwd=tlwd+1 ,  lty=tlty, col=shadow )
                lines( c( subseg[s,"start.pos"], subseg[s,"end.pos"] ), c( subseg[s,segcolumn], subseg[s,segcolumn] ) , 
                   lwd=tlwd ,  lty=tlty, col=segcol[sel][s] )
              }
        abline( v=c(subseg[s,"start.pos"], subseg[s,"end.pos"]), lty=3 )
        mid<-( subseg[s,"start.pos"] + subseg[s,"end.pos"] )/2
        le <- (  subseg[s,"end.pos"] -subseg[s,"start.pos"] )
        if( any( colnames(segments)== "labels"  ) & le>1000000 ){
          axis( 3, at=mid, labels=subseg[s,"labels"] , las=2, cex.axis=1 )
         }
        }

       }

       if( !is.null(n.rc.seg)  ){
         sel<-n.rc.segments$chrom ==chr
         if( sum(sel)>0 ){ 
          subsegn <- n.rc.seg[ sel,]
           for( s in 1:nrow(subsegn)){
            lines( c( subsegn[s,"start.pos"], subsegn[s,"end.pos"] ), 2* c( subsegn[s,"mean"], subsegn[s,"mean"] ) , 
                  lwd= nlwd  , col=ncol, lty=nlty  )
           }
          }
       }

 
  
    
# place holder for allelic ratio graph
    if( !is.null(ar) ){
     sel<-ar$CHR==chr 
       if( sum(sel)>0 ){
        subar<-ar[sel,] 
        subar$ar<-subar[,3]/(subar[,3]+subar[,4])
        hh<-hist2d ( subar$POS, subar$ar , nbins=c( floor( nrow(subar)/50 ), 20 ) , show=F  )
     
        par( mar=c(3,5,0.5,1), xaxs="i", yaxs="i"    )

        plot.new()
        plot.window( xlim=XLIM, ylim=c(0,1)  )
 
        title( ylab="AR") 

        for( s in 1:nrow(subseg)){
           from<- subseg[s,4]
          to<- subseg[s,5]
          sel<-hh$x>=from & hh$x<=to 
          if( sum(sel)>1) {
            if( sum( hh$counts[sel,] != 0  )> 0 )
               image( hh$x[sel], hh$y,  log( hh$counts[sel,])  , col=gray( 1-seq(0,1,.05)  ) ,add=T )
          }
        }

       axis( 1, at=seq(0,300000000, 10000000), labels=as.character( seq(0,300000000, 10000000)/1000 ) )
       axis(2, at=c(0,.5,1) )
       if( !is.null( annotation ) ){
         selectch<- annotation[,1]==chr
         if( sum(selectch )>0 ){ 
           axis( 1, at=annotation[selectch,2], labels=annotation[selectch,3], las=2, cex.axis=cex.annotation )
         }
       }
       
       
       box()
      } else {
       
        plot.new()
        plot.window( xlim=XLIM, ylim=c(0,1)  )
        axis( 1, at=seq(0,300000000, 10000000), labels=as.character( seq(0,300000000, 10000000)/1000 ) )
        axis(2, at=c(0,.5,1) )
        box()

     }

   }  else {
        plot.new()
        plot.window( xlim=XLIM, ylim=c(0,1)  )
        axis( 1, at=seq(0,300000000, 10000000), labels=as.character( seq(0,300000000, 10000000)/1000 ) )
        axis(2, at=c(0,.5,1) )
        box()
      }

 




}



  if( !is.null(file) )
    dev.off() 

}
