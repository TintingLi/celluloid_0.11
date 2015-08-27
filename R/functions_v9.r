

##########################
# THIS USES the copynumber packages

DEBUG<-FALSE

## when calling survey, I want to select a subset of selectedPeaks only and try to fit 2 subclones only 


#######################################################
# this is relative read count. 
  
# This is where I expect to see copy number peaks, 
# given S ( where S=1/ploidy, ploidy being the relative abundance of DNA, S=2 being normal cell)
# T=c( N, T1,T2,...) are subclone fractions, sum(T)=1, N is % normal cells
# C=c( 2, c1,c2,...) are copy number in each clones.
# parameters are S and T


## I divide by 2 so that S is at the position of two copies in tumour
## this assumes that the normal has 2 copies; normal should be inspected
## and segments not consistent with 2 copies be excluded.
## At this moment there are no code to do this, it should be done in 
## the wig files? 

## The parameter vector should sum to 1, it is of length nsubcl+1 

RCC <- function(C,  S, t ){

 if( abs( sum(t) -1 )>1e-6 ){ stop("parameter vector does not sum to 1 (RCC)") }
 return( S*sum( t*C )/2 )

}


#########################################################
# This function creates a matrix (cn) of all copy number possibilities in all normal and subclones
# It creates a global cn object, and a global cnw object (which are the sum of the largest copy
# number differences in subclones + the smallest copy number difference in subclones.
# It can be used as weights. 

# only those configurations are retained: lines with a maximum copy 
# number of maxc (default 2 ) and
# such that there are no more than maxsubcldiff difference in copy 
# number between any two subclones.

prepCN.old<-function( maxc=2 , nsubcl=2 , maxsubcldiff= 2 , maxlines=10000  ){

 if( maxc^nsubcl > maxlines ){ stop( "maxc^nsubcl too large") }

 cn<- expand.grid( rep( list( 0:maxc ), nsubcl ) )
 names(cn)<-paste( "T", 1:nsubcl , sep="")

 if( nsubcl>1 ){
  if( maxsubcldiff>=1 | maxsubcldiff==0  ){
   sel<- apply( cn,1, function(x){ max( dist( x ) ) } ) <= maxsubcldiff 
   cn<-cn[sel,] 
  } else { 
    # if maxsubcl<1 
    # one subclone can have a ratio of copy number of no more than 1/maxsubcldiff
    # relative the copies of the other; 
    # ( when maxsubcldiff < 1)
    # exception made for 0 and 1 copies in one, that are ignored 
    sel<-apply( cn,1, function(x){ mn<-max(1,min(x)) ; mx<-max(x); return( mn/mx>=maxsubcldiff ) } ) 
    cn<-cn[sel,]
  }
 }

 if( nsubcl > 2 & FALSE  ){
   # NOT IMPLEMENTED OR TESTED YET
   # only one subcl can differ from all others (exept when copy number is less than 3 in all?) 
   sel <- apply( cn, 1, function( x ){all(x<3 ) | length( unique(x) ) <=2 } )
   cn<- cn[sel,] 
 }

 # In normal (N), all segments are assumed to have two copies 
 cn <- data.frame( N=2, cn )
 
 # this contains the sum of the max difference between subclones + the min difference.
 # this can be used as weight, it's a simple measure of variability
 
 if(nsubcl>1 ){
  cnw<-as.numeric( apply( cn[,2:ncol(cn)], 1, function(x){ dx<-dist(x);return( max(dx)+min(dx) ) } )  )
  cnw<-cnw+1 
 } else {
  cnw<-rep(1, nrow(cn))
 }

# other possibilities, not tested 
# cnw<- 1/cnw
# cnw<-100* cnw/(sum(cnw) ) 

# making these objects global, so accessible from other functions.
  cnw<<-cnw 
  cn<<- cn

  invisible(cn)

}



#################################################################################################

matpat <-function( cnline ){

 # takes a cn lines (number of copies in each normal+subclones) and returns 
 # the number of "maternal/paternal" combinations

 # in normal, assumed to find 1 and only 1 copy of the reference allele
 li<- list( 1 )
 for( i in 2:length(cnline)  ){
   li[[i]]<- 0:as.numeric(cnline[i] )
 }
 tmp <- expand.grid( li ) 
 out<-c()
 for( i in 1:length(cnline) ){
  out<- cbind( out, tmp[,i] )
  out<- cbind( out, as.numeric(cnline[i])-tmp[,i] )
 }
 out<-as.data.frame(out) 
 colnames(out)[seq(1,ncol(out),2)]<- paste( "m", 0:(length(cnline)-1), sep="" )
 colnames(out)[seq(2,ncol(out),2)]<- paste( "p", 0:(length(cnline)-1), sep="" )
 return(out) 
}


##############################################

# cn has already been called and prepped
# this returns  a list of AR corresponding to each cn lines 
# t are the normal/subclone proportions

# the t parameters sum to 1 

prepAR<-function(t, preserveMatPatDiff=T , preserveMaxSubClDiff=T ){
 
 if( abs( sum(t) -1 )>1e-6 ){ stop("parameter vector does not sum to 1 (prepAR)") }

 if(!preserveMatPatDiff & FALSE ){                                                                                                     
# WHY IS THIS HERE?
 apply( cn, 1, function(x){                                                                                                    
   tmp<-eg( x )                                                                                                                
   e<-t %*% t(tmp) 
   e<-e/sum( t*x )
   return(sort(e) ) 
  } ) # END APPLY   

 }

 tmp<- apply( cn, 1, function(x){ 

   mp<-matpat( x )
 
   m<-mp[ seq( 1,ncol(mp),2 ) ]
   p<-mp[ seq( 2,ncol(mp),2 ) ]

   # removing lines for which the number of (say) paternal chromosomes
   # differ by more than maxsubcldiff. 
   if( preserveMaxSubClDiff & ncol(cn)>2  ){
     maxsubcldiff <- max( apply( data.frame(cn[,2:ncol(cn)]), 1, dist )  )
     sel<- apply( data.frame( m[,2:ncol(m)]), 1, dist )>  maxsubcldiff     
     sel<- sel | apply( data.frame( p[,2:ncol(p)]), 1, dist )>  maxsubcldiff  
     mp<-mp[!sel,]
     m<-m[!sel,]
     p<-p[!sel,]
   }
     
   if(preserveMatPatDiff & ncol(cn)>2 ){
    dif<-m-p
    sel<-apply( sign(dif), 1, function(x){ length( unique(x) )<=2 } )
    mp<-mp[sel,] 
   }
 
    ar<- apply( mp, 1, function(x){
                  m<-x[ seq( 1,length(x),2 ) ]
                  p<-x[ seq( 2,length(x),2 ) ]
                  e<- sum( t*m/sum( t*(m+p) ) )
                  return(e)
               })

    return( sort( ar ) )
  } )

 
  return( tmp )
}


eg<-function( cnline ){
 # in normal, assumed to find 1 and only 1 copy of the reference allele
 li<- list( 1 )
 for( i in 2:length(cnline) ){
   li[[i]]<- 0:as.numeric(cnline[i] )
 }
 expand.grid( li ) 
}



###############
# 
# this returns the distance between observed (selected) and expected peaks
# for each obseved peak, find the expected peak with the smallest distance 
# on the x-axis (the read count axis), then for that expected peak, compute distance 
# on the y-axis (the allele ratio axis).
#
# THIS IS NOT A CONTINUOUS FUNCTION; CONSEQUENCES OF THAT WITH RESPECT TO THE 
# OPTIM FUNCTION USED WERE NOT WELL INVESTIAATED
#
# weightbydiff is typically a cnw vector created by prepCN, or other vector of 
# the same length. It gives weights to corresponding lines of cn
#
# weightbyarp indicates if the distance on the y-axis should be weighted based on 
# the number of ar peaks (motivation: the more peaks there are on the y-axis, the 
# easier it is to find a good match. 
#

# t sums to 1 
# ... argument passed to prepAR 

# wd is a weight given to each point in a cluster of points

# weighbydiff, weightbyarp are not in use anymore

eoDist<-function( selectedPeaks  ,  S, t , wd,
                  xonly=F , weightbydiff=F , weightbyarp=F, 
                  penaltymultiplier=1000, penaltymultipliery=0 ,  ...  ){

 # first element of t is the % of normal cells 
 nsubcl <- length(t)-1 

 if( abs( sum(t) -1 )>1e-8 ){ stop("parameter vector does not sum to 1 (eoDist)") }

 # in case % normal cells is 0 
 if( t[1]==0 ){ t[1]<-1e-8; t<-t/sum(t) } 
 
 # if selectedPeaks has a third column, it
 # indicates where I might find X copy in tumour. 
 # The third column is vector consisting of 0s and a single value = X
 
 
 if( DEBUG ){ cat( ncol( selectedPeaks ) ) }
 if( ncol( selectedPeaks )==3 ){
   if( any( selectedPeaks[,3]>0 ) ){
     forced<-selectedPeaks[,3] 
   } else {
     forced<-rep(0,nrow(selectedPeaks))
   }
 } else {
   forced<-rep(0,nrow(selectedPeaks))
 }

 # If I am not forcing peaks, then I am only allowing decreasing frequencies ( t[2]>=t[3]>=... ) 
 # all % must be positive
 if( nsubcl>1 ){
  if(  any( t<0 | t>1) |( any( t[2:(nsubcl)]< t[3:(nsubcl+1)] ) & all(forced==0 )  )  ){ return( 66666+S+t[2] ) }
 } 
 if( any( t<0 | t>1 )   ){return(77777+S+t[2])}


 # expected peak locations, one the x-axis (read count) 
 ep<- unique( apply( cn, 1, RCC, S=S, t=t ) ) 
 zero<-apply( data.frame(cn[,2:ncol(cn)]==0) ,1, all )
 one<-apply( data.frame(cn[,2:ncol(cn)]==1) ,1, all )

 # I want to force a coverage of the last chosen point 
 #if( all( ep < max( selectedPeaks$x ) ) ){ return( 88888+S+t[2] ) }
 #if( all( ep > min( selectedPeaks$x ) + .1*( ep[one][1] - ep[zero][1] ) ) ){ return( 99999+S+t[2] ) }

 arp<-prepAR(t, ... ) 

 # 
 if( !weightbyarp | TRUE ){
   wy<-rep( 1, length(arp) )
 } else {
   wy<-unlist( lapply( arp, length ) )
 }
 
 # wy<-wy/sum(wy)
 
 # put weight on the x-axis
 if( !weightbydiff | TRUE  ){
   wx<-rep( 1, length(cnw)  )
 } else {
   wx<- cnw
 }

#  wx<-wx/sum(wx)*length( cnw ) 

 mnx<- apply( data.frame( selectedPeaks[,1] ), 1, function(x){
               min( wx*abs( sweep( data.frame(ep) , MARGIN=2, as.numeric(x), FUN="-") ) )
          } )


 # extrating the ep that are closest to the selectedPeaks
 wh<- apply( data.frame( selectedPeaks[,1] , mnx ), 1, 
             function(x){              
                 which( wx*abs( 
                         sweep( data.frame(ep) , MARGIN=2, 
                                as.numeric(x[1]), FUN="-") )  == x[2] )[1]  } )


 # for these, calculating distances on the y-axis
 mny<-apply( data.frame( selectedPeaks[,1:2], wh ), 1, 
             function(x){ 
                min( abs(  sweep( data.frame( arp[ as.integer(x[3]) ][[1]] ) , 
                                    MARGIN=2, as.numeric(x[2]), FUN="-" )  )  )*wy[ x[3] ]  })


 penalty<-0 
 # a way to force a peak to correspond to X copies is to impose a penalty
 # X is a number like 246, which indicates 2 copies in first subcl, 4 in second and 6 in third 
if( DEBUG ){ cat( forced )  }
if(any( forced > 0 ) ){

  if( sum( forced>0 ) >1 ){ stop("can not force more than one peak") }

  fvalues  <- forced[ forced>0 ] 
 
 # if( any( fvalues <10^(nsubcl-1) ) ){stop("problem with forced peak value")}
  # I am recoding the cn into the 246 system
  recode<-apply( data.frame( cn[, 2:ncol(cn)] ), 1, function( x ){ sum( x*10^((length(x)-1):0) ) } ) 

  if( sum( recode==fvalues )==0 ){
    stop( "value found in third column of selectedPeaks does not correspond to any combination in cn")
  }
 
  fcopy<-c()

  for( f in fvalues ){
   fcopy<-c( fcopy, which( recode == f ) )
  }

  penalty<- penaltymultiplier * sum( abs( selectedPeaks$x[ forced>0 ] - ep[ fcopy ] )  )
  penalty<- penalty + penaltymultiplier * penaltymultipliery * sum(mny[  forced>0  ])


}


 if( !xonly ){
   mn<- sqrt( ( ( mnx )^2 + (mny)^2 )  ) 
 } else {
   mn<- mnx 
 }


 return( sum( wd* mn ) + penalty  )

}


####################################################


# N subclones
# the vector v is c( S, T ), T=C(N,T1,T2,..,T(N-1) )
# v is length 2+nsubcl
# the % of the last subclone is set to 1-sum( v[3:lenght(v)] )
# 
# optionally we can set S to a specific value, in which case (is !is.null(S) )
# the vector v is c(T)
# ... are arguments for eoDist
#
# this is the function to be optimized 
#

# v contains S and "t without its last element"
peakProximity <-function( v, selectedPeaks , wd, verbose=T, npeaks =NULL, Sn=NULL, xonly , ... ){

  if( is.null(Sn) ){
   S<-v[1] 
   t<- tail(v, -1 ) 
   t<- c( t, 1-sum(t) )
  } else {
   t<-v
   if( t[1]<0.001 ){ t[1]<-0.001 }
   t<-c( t, 1-sum(t) )
   S<-Sn/t[1]
  }
  
  totDist <- eoDist( selectedPeaks, S=S, t=t, wd=wd ,xonly=xonly,  ... ) 

  if( is.null( npeaks ) ){ npeaks<- nrow(selectedPeaks)}
  # rescaling because a subset of peaks were used. 
  totDist <- totDist*npeaks/nrow( selectedPeaks )

  if( is.na(totDist) ){ totDist <- 55555 }
 
  if( totDist<  thisMn ){
    thisMn <<- totDist 
  }

  # output 
   ta<- c( totDist, S , t )  
  if( totDist< 99999  ){ 
    tmpglobal<-paramSpace 
    paramSpace <<-rbind( tmpglobal, ta ) 
  }

if( DEBUG ){
  cat(verbose,"\n")
  cat(totDist,"\n")
  cat(thisMn,"\n")
}

  if(  (verbose &  totDist<=thisMn & totDist<99999)  ){ 
   # cat( as.numeric(ta), paste( rownames(selectedPeaks), collapse=",") , "\n", sep="\t"   )

     cat(paste( paste("totDist=",totDist, sep="" ) , paste("S=",S, sep=""), paste( "t=c(", paste( t, collapse=", "), ")" ), sep="; "),"\n")
 

#      cat( as.numeric(ta)  , "\n", sep="\t"   )
  } 
  return( totDist )

}


#############

findEquiPeaks<-function( selectedPoints , minpoints=4  ){
 x<-selectedPoints$x
  
 subset<-expand.grid( rep( list( c(TRUE,FALSE) ), length(x) ) )
 subset<-as.matrix(subset[ apply( subset,1,sum)>=minpoints ,] )

 d<-apply( subset, 1, 
  function(sub){
  xx<-x[sub]
  m<-median( xx[2:length(xx)]-xx[1:(length(xx)-1)]  )
  d<-(xx-xx[1])/m
  return( sqrt(  sum(  ( d-floor( d+.5) )^2 )/sum(sub) )  )
 })


 
 
 return( list(subset= subset[ order(d),], d=sort(d)  ) ) 

}






        





 
#################################################################################################


##############


# generates random starting values, for use with optim()
startOptim <-function( Sfrom, Sto , lowerF, upperF , Sn=NULL ){ 
      nsubcl <- length( lowerF ) 
      thisMn<<-99999; 
      if( !is.null( Sfrom ) ){
        S<-runif( 1,Sfrom, Sto );
      } 
      n<-runif( 1, lowerF[1], upperF[1] );
      if( nsubcl== 1 ){
        t<-1-n
      } else {
        t<-rep(NA,nsubcl)
        for( i in 1:(nsubcl-1) ){
          if( upperF[i+1]-n  < lowerF[i+1] ){stop("startOptim: problem with lowerF boundaries")}
          t[i]<-runif( 1, lowerF[i+1], upperF[i+1]-n )
        }
        t[nsubcl]<- 1-n-sum(t, na.rm=T )
      }
      
      v<-c(n,t);
      v<-v/sum(v)
      if( is.null(Sn) ){
        return( c(S, v[1:(nsubcl)] ))
      } else {
        return( v[1:nsubcl] )
      }
 }


###########

# li is a list that contains output from the optimization function used. 
# this function extracts parameters from each list components
# ploidy is a set of bounds for the tumour ploidy 
getLocalSolutions <-function(li, max=TRUE, pruneS=0.05 , pruneT=0.025 , ploidy =NULL ){

 allParamSpace<<-c()
 if( DEBUG) {print (li) }
  out<-data.frame()
  for( i in 1:length(li)){ 
    tmp<-data.frame( t( c( li[[i]]$value, li[[i]]$par, li[[i]]$subset ) ) )
    names( tmp )<-paste( "X",1:ncol(tmp) , sep="")
    out<-rbind( out, tmp  )
    tmpglobal<-allParamSpace
    allParamSpace<<-rbind( tmpglobal , li[[i]]$paramSpace )
  }
  for( i in 1:(ncol( out )-1 )){
    out[,i]<- as.numeric( as.character(out[,i]) ) 
  }
 
  if( max ){ 
     out<-out[ order(out[,1], decreasing = T ),] 
  } else {
    out<-out[ order(out[,1], decreasing = F),]
  }
               
  if( ncol(out)==4 ){
    out<-cbind( out[,1:3], 1- as.numeric(as.character(out[,3])), out[,4] )
    names(out)<-c("value","S", "N", "T1", "subset")
  } else {
    out<-cbind( out[,1:(ncol(out)-1)], 
                1- apply( out[,3:(ncol(out)-1) ], 1, function(x){sum(as.numeric(x))} ), out[,ncol(out)] )  
    names(out)<-c("value","S", "N", paste( "T", 1:(ncol(out)-4 ), sep="") , "subset")
  }
  if( !is.null(pruneS) ){
    dS<-as.matrix( dist(out[, 2] ))
    dT<-as.matrix( dist(out[, 3:(ncol(out)-2)] ))
     rem<-rep( FALSE, nrow(dS))
     for( i in 1:(nrow(dS)-1) ){
       rem <- rem | c( rep(FALSE, i ) ,  dS[i, (i+1):ncol(dS) ] < pruneS  & dT[i, (i+1):ncol(dT) ] < pruneT )
     }
     out<-out[!rem,]
  }
 if( !is.null(ploidy) ){
   rem<-rep( FALSE, nrow(out))
   for( i in 1:nrow(out) ){
     tploi <- (2/out[i,"S"] - 2*out[i,"N"])/(1-out[i,"N"] )
     cat( tploi, ploidy, "\n")
     rem[i]<- tploi< ploidy[1] | tploi > ploidy[2]
   }
   out<-out[!rem,]
 }
  if( max )
    return(out[ order(out[,1], decreasing=T ), 1:(ncol(out)-1) ] )
  else return(out[ order(out[,1]), 1:(ncol(out)-1) ] )

}

###############

getLocalMinsFromParamSpace<-function( pruneS=0.05 , pruneT=0.05 ,pruneDist=2 ){

  out<-paramSpace[ paramSpace[,1]<1000,] 
  out<-out[ order(out[,1]),]
  out<-out[ out[,1]<pruneDist*min(out[,1]),] 
  # dS<-as.matrix( dist(out[, 2] ))
  dT<-as.matrix( dist(out[, 3:(ncol(out))] ))
  
  i<-1
  done<-FALSE
  while( !done ){
    trim <-c( rep( F,i ) ,  dT[(i+1):(nrow(out)), i]<pruneT    )
    # rbind in case only one line
    out<-rbind( out[!trim,] )
    #dS<-dS[ !trim, !trim]
    dT<-dT[ !trim, !trim] 
    i<-i+1
    if( i >= nrow( out ) ){ done<-TRUE }
  }
  
return( out )
}




###########
# won't document this
forcePeak<-function( sp, at, cp ){
  sp$w<--1
for( i in 1:length(at) ){
  d<- as.matrix( dist( c( at[i] ,sp[,1] ) ) )[,1]
  wh<-which( d[-1]==min(d[-1]) )[1]
  sp$w[wh]<- cp[i]
}
  return(sp)
}


#######################################

# won't document this
showParamSpace<-function( outputlist ){

   plot( paramSpace[,1], paramSpace[,ncol(paramSpace)-2 ], pch='.' , ylim=c(0,2), xlab="S", 
          ylab="dist") 
   abline( h=(mn<-min(paramSpace[,ncol(paramSpace)-2 ]) ) , lty=3 )
   for( i in 1:length(outputlist) ){
    abline( v= outputlist[[i]]$par[3] , lty=3 )
    text( outputlist[[i]]$par[3], outputlist[[i]]$value , i )
   }

}

#######################################

# exepected peak positions, given parameters 
ePeakPos <- function( par=NULL, S=NULL, t=NULL, cn, preserveMatPatDiff=T , preserveMaxSubClDiff=T ){

 if( !is.null(par) ){
   if( !is.null( S ) | !is.null(t) ){
    stop("only need one of par or S/t")
   }
   S<-par[1]
   t<-par[2:length(par)]
   t<-c(t, 1-sum(t))
 } else {
   if( is.null(S) | is.null(t) ){
     stop("need either par or S/t")
   }
 }

 mp<- apply( cn, 1 , matpat )
 #require( data.table) 
 # this returns the number of 
 mp <-as.data.frame(rbindlist(mp))
 m<-mp[ seq( 1,ncol(mp),2 ) ]
 p<-mp[ seq( 2,ncol(mp),2 ) ]


 if(  ncol(cn)>2  & preserveMaxSubClDiff  ){
     maxsubcldiff <- max( apply( data.frame(cn[,2:ncol(cn)]), 1, dist )  )
     sel<- apply( data.frame( m[,2:ncol(m)]), 1, dist )>  maxsubcldiff     
     sel<- sel | apply( data.frame( p[,2:ncol(p)]), 1, dist )>  maxsubcldiff  
     mp<-mp[!sel,]
     m<-m[!sel,]
     p<-p[!sel,]
   }

 if( ncol(cn)>2  & preserveMatPatDiff ){
   dif<-m-p
   sel<-apply( sign(dif), 1, function(x){ length( unique(x) )<=2 } )
   mp<-mp[sel,] 
 }


 # this calculates the allelic ratio for each combinations
 # assuming that the "m" chromosomes carry the ref allele
 
 ar<- apply( mp, 1, function(x){
                  m<-x[ seq( 1,length(x),2 ) ]
                  p<-x[ seq( 2,length(x),2 ) ]
                  e<- sum( t*m/sum( t*(m+p) ) )
                  return(e)
               })

 x<-apply( mp, 1, function(x){ 
                     m<-x[ seq( 1,length(x),2 ) ]
                     p<-x[ seq( 2,length(x),2 ) ]
                     return( sum( t*(m+p) ) )
                  }) 

 x<-x*S/2 

 return( data.frame( mp, x, ar ) )

}


################ 


#################################################################################################

# won't document this 
getSegments<-function( copyAr , tc , gamma=500, kmin=100 ){

    stdchr <- c( paste( "chr",1:22, sep="") , "chrX", "chrY" )

    tmp<- data.frame( chrom= tc$space, pos=start( tc)  , rc= tc$icopy  ) 
    tmp<- tmp[ is.element( tmp$chrom, stdchr ),] 

    tmp$chrom<-factor( tmp$chrom, levels= stdchr  )

    tmp$chrom<-as.numeric( tmp$chrom )

    tmp<-tmp[ !is.na( tmp$rc ),] 
    tmp<-tmp[order( tmp$chrom),] 
 
    tmp.win<-winsorize( tmp )
    iseg<- pcf( data=tmp.win, gamma=gamma , kmin=kmin   , digits=4 )

    iseg$chrom<-paste( "chr", iseg$chrom, sep="" )
    iseg$chrom[ iseg$chrom=="chr23" ]<- "chrX"
    iseg$chrom[ iseg$chrom=="chr24" ]<- "chrY"

    copyAr$iseg <<- iseg 
    invisible(iseg)

}






