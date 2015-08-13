
peakProximity.percent <-function( v, segments , verbose=T, Sn=NULL , ... ){
  
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
  
  totDist <- eoDist( segments , S=S, t=t , ... ) 
  
  if( totDist <  thisMn ){
    thisMn <<- totDist 
  }
  
  # output 
  ta<- c( totDist, S , t )  
  if( totDist < 1  ){ 
    tmpglobal<-paramSpace 
    paramSpace <<-rbind( tmpglobal, ta ) 
  }
  
  if( DEBUG ){
    cat(verbose,"\n")
    cat(totDist,"\n")
    cat(thisMn,"\n")
  }
  
  if(  (verbose &  totDist<=thisMn & totDist<1)  ){     
    cat(paste( paste("totDist=",totDist, sep="" ) , paste("S=",S, sep=""), 
               paste( "t=c(", paste( t, collapse=", "), ")" ), sep="; "),"\n")
  } 
  return( totDist )
  
}


##################################################################
# here segments is an obsect such as t.ar.seg. It must have columns named mean, p and size
# this only works for one clone models 

eoDist.percent<-function( segments ,  S, t , quadratic=TRUE ,...  ){
  
  if( sum( is.element( c("mean","p","size"), colnames(segments) ) )!=3){
    stop("eoDist.percent: segments must have columns named mean, p and size")
  }
  
  seg<-segments[ !apply( is.na(segments[,c("mean","p") ] ), 1, any ), ]
  
  # first element of t is the % of normal cells 
  nsubcl <- length(t)-1 
  
  if( nsubcl > 1 ){ stop("eoDist.percent only implemented for one clone models") }
  
  if( abs( sum(t) -1 )>1e-8 ){ stop("parameter vector does not sum to 1 (eoDist.percent)") }
  
  # in case % normal cells is 0 
  if( t[1]==0 ){ t[1]<-1e-8; t<-t/sum(t) } 
  
  # I am only allowing decreasing frequencies ( t[2]>=t[3]>=... ) 
  # all % must be positive
  if( nsubcl>1 ){
    if(  any( t<0 | t>1) |( any( t[2:(nsubcl)]< t[3:(nsubcl+1)] )   )  ){ return( 1 ) }
  } 
  
  if( any( t<0 | t>1 )   ){return(1)}
    
  epp<-ePeakPos( S=S, t=t, cn=cn, ... )
  
  d<-dist( epp$x )
  xdist <-min( d[d>0] )/2
  
  nseg<-nrow(seg)
  nepp<-nrow(epp)
  
  eoseg<-rbind(as.matrix( seg[,c("mean","p")]  ), as.matrix( epp[,c("x","ar")] ) )
  
  d<-as.matrix(dist(eoseg) )
  d<-d[1:nseg,(nseg+1):(nepp+nseg)] 
  
  mnd<-apply( d, 1, min )
  
  if( quadratic ){
    mnd<-weigth.quadratic( mnd, xdist=xdist)
  } else {
    mnd<-weigth.linear( mnd, xdist=xdist)
  }
  
  glength<-sum( seg$size )
  percent<-sum( mnd * seg$size )/glength
  # so that the function can be minimized 
  return( 1 - percent )
  
}


###########################################################################################################


coverParamSpace.percent <- function(  segments, verbose=T , addToParamSpace=F , control=NULL , 
                             Sfrom=NULL, Sto=NULL, Sn=NULL, 
                             maxc=NULL, maxsubcldiff=NULL , optimFct=2 , lowerF, upperF , 
                             nrep=1 , method=NULL, percentObj=TRUE, ...  ){
  
  objectiveFct<-peakProximity
  if( percentObj ){   objectiveFct<-peakProximity.percent }
  
  if( is.null(method) ){method="L-BFGS-B"}
  
  totDist<<-Inf  
  rownames( segments )<- 1:nrow(segments)
  
  # holds the output of optimization functions
  outputlist <-list() 
  
  nullmaxsubcldiff<-FALSE  
  if( is.null( maxsubcldiff) ){ nullmaxsubcldiff<-TRUE } 
    
  if( is.null(Sn) ){
    if( is.null(Sfrom) | is.null(Sto) ){
      stop("Sfrom, Sto: undefined")
    }
  } else {
    Sfrom<-NULL
    Sto<-NULL
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
    
    if( !exists( "paramSpace"  ) | !addToParamSpace  | nrep>1  ){
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
                     segments=segments , ... ) 
      } else {
        op <- GenSA( par=NULL ,fn=objectiveFct,lower=c( lowerF )  ,   
                     upper=c( upperF ),  control=control , Sn=Sn, 
                     segments=segments , ... ) 
      }
      
      
      if( !is.null( Sn ) ){ op$par<- c( Sn/op$par[1], op$par ) }
      
      outputlist[[currentrep]] <- op
      outputlist[[currentrep]]$start<-NULL
      outputlist[[currentrep]]$subset<- paste( subset, collapse=",")
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      outputlist[[currentrep]]$paramSpace<-paramSpace
    }
    
    ###########################################
    # use optim, a number of times; not grid 
    if( length(optimFct)==1 & optimFct[1]==1 ){ # use optim
      if( is.null(control) ){ 
        stop("Error: must specify control. See ?optim.") 
      }
      # starting values  
      start<-startOptim( Sfrom, Sto, nsubcl , forced=FALSE , Sn=Sn )
      # no need to split, as parameter list is defined in startOptim
      op<- optim( par=start , fn=objectiveFct, segments=segments[subset,]  ,
                  verbose= T , control=control, npeaks=nrow( segments), Sn=Sn ,... ) 
      
      if( !is.null( Sn ) ){ op$par<- c( Sn/op$par[1], op$par ) }
      
      outputlist[[currentrep]]<-op
      outputlist[[currentrep]]$start<-start
      outputlist[[currentrep]]$subset<- paste( subset, collapse=",")
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      outputlist[[currentrep]]$paramSpace<-paramSpace
      
    }
    
  } # for( currentrep in 1:nrep)
  
  # grid search, grid consist of starting values
  # if optimFct > 2 then take that many grid points for each parameters
  # does not need nrep to be defined, the number of rep is determined by number of starting points
  
  if( (length(optimFct)==1 & optimFct[1]>2 ) | length(optimFct)>1 ){
    if( length(optimFct)>1 ){
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
          op<- optim( par=start , fn=objectiveFct, segments=segments,  
                      lower=c( Sfrom , lowerF )  ,   upper=c( Sto, upperF ),
                      method=method, verbose= T , control=control, npeaks=nrow( segments) , ... ) 
        } else {
          op<- optim( par=start , fn=objectiveFct, segments=segments[subset,],  
                      lower=c( lowerF )  ,   upper=c( upperF ),
                      method=method, verbose= T , control=control, npeaks=nrow( segments) , Sn=Sn,...)
        }
        if( !is.null( Sn ) ){ op$par<- c( Sn/op$par[1], op$par ) }
        
        outputlist[[currentrep]]<-op
        outputlist[[currentrep]]$start<-start
        outputlist[[currentrep]]$subset<- paste( subset, collapse=",")
        outputlist[[currentrep]]$maxc<-maxc
        outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
        outputlist[[currentrep]]$paramSpace<-paramSpace
        outputlist[[currentrep]]$cn<-cn
      }
    }
  }
  
  
  invisible(outputlist) 
  
}


################################


weight.linear <-function(di, xdist){ 
  we <- 1-di/xdist;
  we[ di>xdist ]<-0
  return(we)
}

weight.quadratic <-function(di, xdist ){ 
  b<- -2/xdist 
  a<- -b/(2*xdist)
  c<- 1
  we<- a*di^2+b*di+c  
  we[ di > xdist ] <- 0 
  return(we)
}