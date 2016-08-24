
showSTumourProfile<-function(scopyAr, nx=200, ny=50, maxPoints=50000, selected=NULL, 
                            flatten=.25, xlim=c(0,2), nlev=50 , noise=NULL  ){
 #require(MASS)
 if( is.finite(maxPoints) ){
   sa<-sample( 1:nrow(scopyAr), min( maxPoints, nrow(scopyAr)) , replace=F )
 } else { 
   sa<-1:nrow( scopyAr )
 }
 
 if( !is.null( noise ) ){
   noise<- rnorm( length( sa ) , 0, noise) 
 } else { noise<-rep(0, length(sa) ) } 
 
 cntr <-kde2d( scopyAr[sa,"copy"] +noise , scopyAr[sa,"sar" ]  , n=c( nx, ny ) , lims=c( xlim, 0,1 ) )
 plot(  scopyAr[sa,1] + noise , scopyAr[sa,2] , pch='.', col="gray", ylim=c(0,1), 
         xlim=xlim , xaxt = "n" , xlab="Copy number" , ylab="SV AR (%AltAlleles)" ) 
 image(  cntr$x, cntr$y, cntr$z^flatten ,add=T , col=terrain.colors(50) )
 contour( cntr$x, cntr$y, cntr$z^flatten , xlim=xlim , nlev=nlev  , add=T )

}


###################

############# NEED TO CHECK WHICH ONE IN SAR IS THE REF ALLELE

# copySAr will have the same columns as copyAr
prepSomaticCopyAr<-function( tumourrangedata, sar ){

  tmp<-tumourrangedata[1,]
  binlength<-end( tmp )-start(tmp) + 1
  
  ars <- data.frame( chpo=paste( tumourrangedata$space, start(tumourrangedata), sep="-") ,copy= tumourrangedata$smcopy   ) 
  alratio<-sar
  chr<- alratio[,1]
  pos<- alratio[,2]
  alratio$chpo <- paste( chr , floor( pos/binlength )*binlength+1  , sep="-") 
  #alratio$nar <- alratio[,5]/( alratio[,5]+alratio[,6])
  alratio$sar <- alratio[,4]/( alratio[,3]+alratio[,4])
  ars <-merge( ars , alratio, by="chpo") 
  
  ars<-ars[ sel<- !apply( is.na(ars), 1, any ), ] 
  return( data.frame( ars[, c("copy","sar")], chr=ars[,3], pos=ars[,4], ref=ars[,5], alt=ars[,6])  ) 
}



###################################################################

# x is selectedPoints 

addSomaticLabels_ALT<-function(x , eSP,  ep , all=F ){
  if( !all ){
    for( i in 1:nrow(x) ){
      # select the point from ep that is closest to x 
      names(x)<-c("x","y")
      d<-as.matrix( dist( rbind(x[i,], ep[,1:2] ) )  )
      sel<-which( d[,1][-1] == min( d[,1][-1] ) )[1]
      
      #pick in eSP the positions that match m1,p1,m2,p2 in ep
      cn<-ep[sel,c("m1","p1","m2","p2")]
      sel<- apply( apply( eSP[, c("m1","p1","m2","p2") ] , 1, function(x){x==cn} ) , 2, all)
      
      sub<-eSP[sel,] 
      
      for( i in 1:nrow( sub ) ){
        chrlab<- paste("",paste( paste( sub[i ,c("m1","p1") ] , collapse="" ),   
                                 paste( sub[i ,c("m2","p2") ], collapse="" ),sep="/"), sep="")
        somlab<- paste("",paste( paste( sub[i ,c("sm1","sp1") ], collapse="" ),   
                                 paste( sub[i ,c("sm2","sp2") ], collapse="" ),sep="/"),sep="")
        legend( sub[i,"x"],  sub[i,"y"], 
                legend=  c( chrlab, somlab ),
                cex=.5, bg="white" )
      }
    }
  }
}


############

addSomaticLabels<-function(  eSP , manual=T  ){
  
  nc<-ncol(eSP)
  nsubcl <- length( grep( "sm", names( eSP ) ) )
  if( manual ){ 
    
    id<-identify( eSP$x, eSP$sar , n=1 , labels=""  )
    while( length(id) > 0 ){
      # THIS WON"T WORK FOR SUBCLONES
      m<-eSP[id,][seq(3, 2*nsubcl+1  ,2) ]
      p<-eSP[id,][seq(4, 2*nsubcl+2  ,2) ]
      
      chrlab <-  paste( m+p , collapse="/") 
      somlab<- paste( eSP[id, grep( "sm", names(eSP) )  ] , collapse="/" )
      
      legend( eSP[id,"x"],  eSP[id,"sar"], 
              legend=  c( chrlab, somlab ),
              cex=.5, bg="white" )
      
      id<-identify( eSP$x, eSP$sar , n=1 , labels=""  )
      
    }
    
    
  } else {
    for( id in 1:nrow( eSP ) ){
      m<-eSP[id,][seq(3, 2*nsubcl+1  ,2) ]
      p<-eSP[id,][seq(4, 2*nsubcl+2  ,2) ]
      
      chrlab <-  paste( m+p , collapse="/") 
      somlab<- paste( eSP[id, grep( "sm", names(eSP) )  ] , collapse="/" )
      
      legend( eSP[id,"x"],  eSP[id,"sar"], 
              legend=  c( chrlab, somlab ),
              cex=.5, bg="white" )
    }
  }
}





plotSomaticModelPeaks<-function( S, t , cn, col="red",cex=1,lwd=3 , plot=T , eSP=NULL ){
  
  
  epp <- ePeakPos(S=S, t=t ,cn=cn )

  ep<- ( apply( cn, 1, RCC, S=S, t=t ) ) 
  d<-data.frame(cn,ep)
  d<-d[ order( d$ep),] 
  if( plot ){
    s<-seq(1,length(ep),2)
    axis( 1, at=d$ep[s] , labels=apply( as.matrix(d[,2:ncol(cn)]), 1, paste, collapse =".")[s] , las=2  , cex.axis=.75) 
    s<-seq(2,length(ep),2)
    axis( 3, at=d$ep[s] , labels=apply( as.matrix(d[,2:ncol(cn)]), 1, paste, collapse =".")[s] , las=2  , cex.axis=.75  ) 
  }
  sel<- apply( epp, 1, function(x){  return( all( x[seq( 3,length(x)-2, 2 ) ]==x[3] &  x[seq( 4,length(x)-2, 2 ) ]==x[4]  ) ) } ) 
  
  x<-unique( epp[sel,"x"] )
  abline( v=x, col=col)

  if( is.null( eSP) )
    eSP<-eSomaticPeakPos( t,epp  )
  
  if( plot ){
    xy<- unique( eSP[,c("x","sar")  ] ) 
    points( xy[,1], xy[,2] , pch=21,  col=col, cex=cex, lwd=lwd  )
  }
  
  invisible(eSP) 
  
}




eSomaticPeakPos <-function( t , epp ){
  
  
  nsubcl <- (ncol(epp)-4 )/2 
  output<-  apply( epp , 1, 
                   function(x){
                     x<-as.numeric(x)
                     mat<-x[seq(3,ncol(epp)-2 ,2)]
                     pat<-x[seq(4,ncol(epp)-2 ,2)]
                     copy<- mat+pat 
                     li<- list()
                     for( i in 1:nsubcl  ){
                       li[[i]]<- 0:copy[i]
                     }
                     grid <- expand.grid( li ) 
                     sar<- apply( grid, 1, 
                                    function(g){ 
                                      sum( g*t[2:(nsubcl+1)] )/(2*t[1]+sum( t[2:(nsubcl+1)]*copy ) )
                                    }
                     )
                     d<-data.frame( matrix( x, ncol=ncol(epp), nrow=nrow(grid), byr=T ) , grid, sar )
                     return(d) 
                   }
  )
  
  # require( data.table) 
  output <-as.data.frame(rbindlist(output))
  names( output )[ 1:ncol(epp) ]<-names(epp)
  names( output )[ (ncol(epp)+1):ncol(output) ] <- c( paste( "sm", 1:nsubcl , sep="") ,"sar") 
  return(output)
  
}














