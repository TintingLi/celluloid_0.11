
  
arInSeg<-function(seg, ar, minhet=50 , maxhet=5000 , control=list( maxit=1000 ),
                          tumourrangedata=NULL, maskmap=.8   ){

 # take a het, has n reads.  x of them are showing the ref allele. 
 # p is the average proportion of chromosome bearing the ref allele (theor speaking,
 # the average over all cells). 
 # Modeling x after a poisson. binomial doesnt work. the higher variance of poission 
 # makes it more appropriate
 # q is nuisance paramters, it estimate the proportion of hets having ref allele falling
 # on the same parental chromosome.

  ar.tmp<-ar

  if( maskmap>0 & !is.null( tumourrangedata ) ){
    binlength<- end( tumourrangedata[1,] ) -start( tumourrangedata[1,] ) +1 
    # in which bin does the het belong to 
    bin<-paste( ar[,1], binlength*floor( ar[,2]/binlength )+1, sep="-") 
    sel<- tumourrangedata$map < maskmap 
    mask <- paste( tumourrangedata$space, start(tumourrangedata), sep="-")[sel]
    ar.tmp<-ar[!is.element( bin, mask ),] 
  }
 

dz<-function( z, C, p ){
 return( dpois( z, C*p ) + dpois( z,C*(1-p) ) - ppois( z, C*p )*ppois( z, C*(1-p) ) +  ppois( z-1, C*p )*ppois( z-1, C*(1-p) ))
}


zlik<-function( pr, z , C ){
#  p<-exp(pr)/(1+exp(pr))
  p<-pr
  return(   -sum(  log( dz( z, C , p ) ) )  )
}

  
 tmpseg<-seg
 tmpseg$p <- NA

 for( i in 1:nrow(tmpseg) ){cat("#########",i,"\n")
  chr<-as.character(tmpseg$chrom[i])
  from<-tmpseg$start.pos[i]
  to<-tmpseg$end.pos[i]
  sel<-ar.tmp[,1]==chr & ar.tmp[,2] >=from & ar.tmp[,2] <=to 
  if( sum(sel)>= minhet ){
    x<-ar.tmp[sel, 3]
    y<-ar.tmp[sel, 4]
    bpx<-boxplot(x, plot=F)
    bpy<-boxplot(y, plot=F)
    rem<-is.element( x, bpx$out ) | is.element( y, bpy$out ) | x+y==0 
    x<-x[!rem] ;y<-y[!rem]; n<-x+y    
    sel<-sample( 1:length(x),min( length(x), maxhet ) , replace=F ) 
    x<-x[sel]; n<-n[sel]
    x<-apply( data.frame( x,n-x), 1, min )
    pr<- max( 0.005,min( tmpseg$meanar[i], 1-tmpseg$meanar[i])  )
    cat(i," ",pr,"\n")
    bestmo <-tryCatch(     stats::optimize( f=zlik, interval=c(0,.5), z=x, C=n ) , error=function(e){})
    if( !is.null(bestmo) ){
     print((bestmo))
     est<-bestmo$minimum
     tmpseg$p[i]<- est 
   }
   }
    print( tmpseg[i,] ) 
#if(0){
#    plot( tmpseg[1:i,c("meanar","p")], xlim=c(0,1), ylim=c(0,1) )
#    abline(0,1); abline(1,-1)
#    abline( h=.5); abline( v=.5)
#}
 }


  return(tmpseg )

}



