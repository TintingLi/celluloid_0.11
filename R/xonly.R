# this function is 
# only used if nsubcl=1 no need to do a check unless I modify functions below 

ePeakPos.xonly <- function( par=NULL, S=NULL, t=NULL, cn ){
  
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
  
  mp<- apply( cn, 1 , celluloid:::matpat )
  #require( data.table) 
  # this returns the number of 
  mp <-as.data.frame(rbindlist(mp))
  
  # this calculates the allelic ratio for each combinations
  # assuming that the "m" chromosomes carry the ref allele
  
  # faking ar to be .5 
  ar<- rep( 0.5, nrow( mp ) ) 
  
  x<-apply( mp, 1, function(x){ 
    m<-x[ seq( 1,length(x),2 ) ]
    p<-x[ seq( 2,length(x),2 ) ]
    return( sum( t*(m+p) ) )
  }) 
  
  x<-x*S/2 
  
  return( data.frame( mp, x, ar ) )
  
}

#########################################
#########################################
#########################################
#########################################

peakProximity.xonly  <-function( v, copyArTable  , verbose=T, Sn=NULL , ... ){
  
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
  
  totDist <- eoDist.xonly( copyArTable  , S=S, t=t , ... ) 
  
  if( totDist <  thisMn ){
    thisMn <<- totDist 
  }
  
  # IN ALL OUTPUT I USE PERCENT CAPTURE, BUT THE CODE 
  # MINIMIZED 1-PERCENTCAPTURE
  # output 
  ta<- c( 1-totDist, S , t )  
  if( totDist < 1  ){ 
    tmpglobal<-paramSpace 
    paramSpace <<-rbind( tmpglobal, ta ) 
  }
  
  if( DEBUG ){
    cat(verbose,"\n")
    cat(totDist,"\n")
    cat(thisMn,"\n")
  }
  
  if(  (verbose &  totDist<=thisMn )  ){     
    cat(paste( paste("value=",1-totDist, sep="" ) , paste("S=",S, sep=""), 
               paste( "t=c(", paste( t, collapse=", "), ")" ), sep="; "),"\n")
  } 
  return( totDist )
  
}

############################
############################
############################
############################

##################################################################

eoDist.xonly <-function( copyArTable  ,  S, t , quadratic=TRUE ,...  ){
  
  # first element of t is the % of normal cells 
  nsubcl <- length(t)-1 
  
  if( nsubcl > 1 & FALSE  ){ stop("eoDist.xonly only implemented for one clone models") }
  
  if( abs( sum(t) -1 )>1e-8 ){ stop("parameter vector does not sum to 1 (eoDist.xonly)") }
  
  # in case % normal cells is 0 
  if( t[1]==0 ){ t[1]<-1e-4; t<-t/sum(t) } 
  if( t[1]==1 ){ t[1]<-.9999; t<-t/sum(t) }
  # I am only allowing decreasing frequencies ( t[2]>=t[3]>=... ) 
  # all % must be positive
  
  if( any( t<=0 ) ){ return( 1 - min(t) ) }
  if( any( t>=1 ) ){ return( max( t) ) }
  if( S> max( copyArTable$mean , na.rm=T ) ){ return( 1+S ) }  
  
  epp<-ePeakPos.xonly( S=S, t=t, cn=cn, ... )
  
  
  
  nseg<-nrow(copyArTable)
  nepp<-nrow(epp)
  
  eoseg<-c( copyArTable$mean, epp[,c("x")] )
  
  sel0 <- apply( epp[, 3:(ncol(epp)-2 )]==0 , 1, all )
  if( ncol(epp)-2 > 4 ){ # if there are subclones
    sel1 <- apply( epp[, seq(3,(ncol(epp)-2 ),2)]==1 , 1, all ) & 
      apply( epp[, seq(4,(ncol(epp)-2 ),2)]==0 , 1, all )
  } else {
    sel1<- epp[,3]==1 & epp[,4]==0 
  }
  # the x distance between 0 copy and 1 copy
  xdist<- (epp[sel1,"x"]-epp[sel0,"x"])/(nsubcl*2)
  
  d<-as.matrix(stats::dist(eoseg) )
  # distances between segments and epp (only on the x-axis)
  d<-d[1:nseg,(nseg+1):(nepp+nseg)] 
  
  # for each segment, take the min one
  mnd<-apply( d, 1, min )
  
  if( quadratic ){
    mnd<-celluloid:::weight.quadratic( mnd, xdist=xdist)
  } else {
    mnd<-celluloid:::weight.linear( mnd, xdist=xdist)
  }
  
  glength<-sum( copyArTable$size )
  percent<-sum( mnd * copyArTable$size )/glength
  
  # so that the function can be minimized 
  return( 1 - percent )
  
}


############################
############################
############################
############################


coverParamSpace.xonly <- function(  copyArTable  , verbose=T , addToParamSpace=F , control=NULL , 
                                    Sfrom=NULL, Sto=NULL , 
                                    maxc=NULL, maxsubcldiff=NULL , optimFct=2 , lowerF, upperF , 
                                    nrep=1 , method=NULL, ...  ){
  Sn<-NULL 
  
  objectiveFct<-peakProximity.xonly 
  
  if( is.null(method) ){method="L-BFGS-B"}
  
  totDist<<-Inf  
  
  
  # holds the output of optimization functions
  outputlist <-list() 
  
  nullmaxsubcldiff<-FALSE  
  if( is.null( maxsubcldiff) ){ nullmaxsubcldiff<-TRUE } 
  
  if( is.null(Sfrom) | is.null(Sto) ){
    stop("Sfrom, Sto: undefined")
  }
  
  k<-0
  
  if( length(optimFct)==1 &  optimFct[1]==2 & !is.null(nrep) ){ 
    if( nrep>1 ){
      cat("GenSA will search for a global minimum, resetting nrep to 1\n")
      nrep<-1
    }
  }
  
  
  if( length(optimFct) > 1 & !is.null(nrep) ){ 
    if( nrep>1 ){
      cat("Grid search, resetting nrep to 1\n")
      nrep<-1
    }
  }
  
  if( length(optimFct)==1 &  optimFct[1]==1  & is.null(nrep) ){ stop("nrep parameter missing\n") }
  
  # if null, defining it for the loop
  if( is.null(nrep) ){ nrep=1 }
  
  
  ########################################################################
  
  # preparing the cn data.frame containing the allowed copy number configuration is normal and all subclones. 
  # if maxc is NULL, I am using the existing one
  if( is.null( maxc )){
    if( !exists("cn") ){ stop("can'r find cn\n") }
    if( ncol( cn )!= length( lowerF )+1 ){stop("cn doesn't have the right number of subclones\n") }
  } else { 
    if( nullmaxsubcldiff ){ 
      cat("setting maxsubcldiff to ", maxc,"\n") 
      maxsubcldiff<-maxc
    } 
    # lowerF parameter defines the number of subclones 
    nsubcl<- length( lowerF ) 
    # the function greates a global cn, so can be re-used and refered to
    prepCN( maxc , nsubcl, maxsubcldiff )
  }
  
  ########################################################################
  
  
  # current minimum distance between expected and observed peaks
  
  for( currentrep in 1:nrep ){
    
    if( !exists( "paramSpace"  ) | !addToParamSpace   ){
      # global, holds parameters in each (or best so far?) iterations. Will be a data.frame
      paramSpace <<- c()
    }    
    thisMn<<-99999
    
    ########################################
    # use GenSA, simulated annealing.It will search for a global min, nrep was forced to 1 above. 
    if( length(optimFct)==1 & optimFct[1]==2 ){
      #require(GenSA)
      if( is.null(control) ){ 
        stop("Error: must specify control. See ?GenSA.") 
        control=list( maxit=500 ) 
      }
      # I need to split because the vector parameter is not the same
      if( is.null( Sn ) ){
        op <- GenSA( par=NULL ,fn=objectiveFct,lower=c( Sfrom , lowerF )  ,   
                     upper=c( Sto, upperF ),  control=control , 
                     copyArTable=copyArTable , ... ) 
      } else {
        op <- GenSA( par=NULL ,fn=objectiveFct,lower=c( lowerF )  ,   
                     upper=c( upperF ),  control=control , Sn=Sn, 
                     copyArTable=copyArTable , ... ) 
      }
      
      
      if( !is.null( Sn ) ){ op$par<- c( Sn/op$par[1], op$par ) }
      # because minimization is done, and focus is on fraction capture:
      op$value<-1-op$value 
      
      outputlist[[currentrep]] <- op
      outputlist[[currentrep]]$start<-NULL
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      outputlist[[currentrep]]$subset<-"ALL"
      outputlist[[currentrep]]$paramSpace<-paramSpace
    }
    
    ###########################################
    # use optim, a number of times; not grid 
    if( length(optimFct)==1 & optimFct[1]==1 ){ # use optim
      if( is.null(control) ){ 
        stop("Error: must specify control. See ?optim.") 
      }
      # starting values  
      start<-celluloid:::startOptim( Sfrom, Sto, lowerF=lowerF, upperF=upperF , Sn=Sn )
      cat( "replicate: ",currentrep,"/",nrep, "; start:" ,start,"\n")
      if( is.null(Sn)){
        op<- optim( par=start , fn=objectiveFct, copyArTable=copyArTable,  
                    lower=c( Sfrom , lowerF )  ,   upper=c( Sto, upperF ),
                    method=method, verbose= T , control=control , ... ) 
      } else {
        op<- optim( par=start , fn=objectiveFct, copyArTable=copyArTable,  
                    lower=c( lowerF )  ,   upper=c( upperF ),
                    method=method, verbose= T , control=control, Sn=Sn,...)
      }
      
      
      
      if( !is.null( Sn ) ){ op$par<- c( Sn/op$par[1], op$par ) }
      # because minimization is done, and focus is on fraction capture:
      op$value<-1-op$value 
      
      outputlist[[currentrep]]<-op
      outputlist[[currentrep]]$start<-start
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      outputlist[[currentrep]]$subset<-"ALL"
      
      outputlist[[currentrep]]$paramSpace<-paramSpace
      
    }
    
  } # for( currentrep in 1:nrep)
  
  # grid search, grid consist of starting values
  # if optimFct > 2 then take that many grid points for each parameters
  # does not need nrep to be defined, the number of rep is determined by number of starting points
  
  if( (length(optimFct)==1 & optimFct[1]>2 ) | length(optimFct)>1 ){
    #if( length(optimFct)>1 ){
    if( length(optimFct)==1 ){ optimFct<-rep( optimFct[1], length(upperF)+1) }  
    if(is.null(Sn)){
      grid<-list( seq( Sfrom, Sto, len=2*optimFct[1]+1 ) )
      grid[[1]]<-grid[[1]][  seq(2,length(grid[[1]]),2)  ]
      for( co in 1:length( lowerF)  ){
        grid[[co+1]]<-seq( lowerF[co], upperF[co], len=2*optimFct[co+1]+1 )
        grid[[co+1]]<-grid[[co+1]][ seq(2,length(grid[[co+1]]),2)  ]
      }
    } else {
      grid<-list()
      for( co in 1:length( lowerF)  ){
        grid[[co]]<-seq( lowerF[co], upperF[co], len=2*optimFct[co]+1 )
        grid[[co]]<-grid[[co]][ seq(2,length(grid[[co]]),2)  ]
      }
    }
    egrid<-expand.grid( grid )
    for( currentrep in 1:nrow(egrid) ){
      # Overwrites these global par for each grid point.
      
      # paramSpace<<-c()
      
      thisMn<<-99999; 
      start<-as.numeric( egrid[currentrep,] )
      cat( start ,"\n") 
      if( is.null(control) ){ 
        stop("Error: must specify control. See ?optim.") 
      }
      if( is.null(Sn)){
        op<- optim( par=start , fn=objectiveFct, copyArTable=copyArTable,  
                    lower=c( Sfrom , lowerF )  ,   upper=c( Sto, upperF ),
                    method=method, verbose= T , control=control , ... ) 
      } else {
        op<- optim( par=start , fn=objectiveFct, copyArTable=copyArTable,  
                    lower=c( lowerF )  ,   upper=c( upperF ),
                    method=method, verbose= T , control=control, Sn=Sn,...)
      }
      if( !is.null( Sn ) ){ op$par<- c( Sn/op$par[1], op$par ) }
      # because minimization is done, and focus is on fraction capture:
      op$value<-1-op$value 
      
      outputlist[[currentrep]]<-op
      outputlist[[currentrep]]$start<-start
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      outputlist[[currentrep]]$subset<-"ALL"
      
      outputlist[[currentrep]]$paramSpace<-paramSpace
      outputlist[[currentrep]]$cn<-cn
    }
    #}
  }
  
  
  invisible(outputlist) 
  
}



