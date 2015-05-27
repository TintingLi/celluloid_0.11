# lowerF (upperF) gives lower (upper) bounds for the vector of parameters c( N, T1,T2,..,T(N-1) )
# optimFct=2 uses the GenSA function.
# optimFct=1 uses optim
# ... is passed to the function to be minimized 

# from help( is.integer )
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

coverParamSpace<- function(  selectedPeaks, verbose=T , addToParamSpace=F , control=NULL , 
                             Sfrom=NULL, Sto=NULL, Sn=NULL, 
                             maxc=NULL, maxsubcldiff=NULL , optimFct=2 , lowerF, upperF , 
                             nrep=NULL  , usesubsets=NULL , xonly=FALSE, 
                             modeat=NULL, weight=NULL, notSeenPenalty=TRUE  , method=NULL, ...  ){
 
  if( is.null(method) ){method="L-BFGS-B"}

  totDist<<-Inf
  
  if( FALSE ){
    WDTHRESHOLD<-0.05
    if( !xonly ){  
      # selectedPeaks may contain multiple points that are close one another
      # I want to downweight them so that they are not collectively "pulling" a solution too much
      # the weight is just based on the number of points with distance less than 0.05
      # (CONSIDER passing wd or treshold as argument instead?)
      tmp<-selectedPeaks[,1:2]
      sel<-tmp[,2]>.5
      tmp[sel,2]<-1-tmp[sel,2]
      d<-as.matrix(dist(tmp[,1:2] ))
      wd<-1/apply( d< WDTHRESHOLD  , 1, sum )
      wd<-sqrt(wd/sum(wd)*nrow(selectedPeaks))
    } else {
      d<-as.matrix(dist(selectedPeaks[,1] ))
      wd<-1/apply( d< THRESHOLD , 1, sum )
      wd<-sqrt(wd/sum(wd)*nrow(selectedPeaks))
    }
  }
  
  # I DECIDED TO IGNORE WEIGHTS
  wd<-rep(1,nrow(selectedPeaks) ) 

  # Only useful if selectedPeaks is based on individual segments
  # modeat is either NULL or take integer value. The most frequent (mode) x value will be forced to correspond to an 
  # absolute copy num
  if( !is.null( modeat ) ){
   
    if( !is.wholenumber( modeat ) | modeat<0  ){stop("Error: modeat should take an integer value")}
    d<- density( selectedPeaks[,1] , weight=weight )
    mo <- d$x[which( d$y==max(d$y) )]
    d<- as.matrix( dist( c( mo,selectedPeaks[,1] ) ) )[,1]
    wh<-which( d[-1]==min(d[-1]) )[1]
    selectedPeaks$w<- -1
    selectedPeaks$w[wh]<-modeat
    print(selectedPeaks)
  }
  
  # Force an x value to represent two copies. Closest point in selectedPeak will be forced to correspond to
  # two values. NUKING THIS 
#   twoat<-NULL
#     if( !is.null( twoat ) ){
#       d<- as.matrix( dist( c( twoat ,selectedPeaks[,1] ) ) )[,1]
#       wh<-which( d[-1]==min(d[-1]) )[1]
#       selectedPeaks$w<- -1
#       selectedPeaks$w[wh]<-  as.numeric(paste(rep(2, length(upperF) ), collapse="") )
#       print(selectedPeaks)
#     }


  rownames( selectedPeaks )<- 1:nrow(selectedPeaks)
  
  # holds the output of optimization functions
  outputlist <-list() 
  
  nullmaxsubcldiff<-FALSE  
  if( is.null( maxsubcldiff) ){ nullmaxsubcldiff<-TRUE } 
  
  # if( nrow( selectedPeaks )<2 ){ stop("Only 1 selected point found") }
  
  if( is.null(Sfrom) ){
    stop("Sfrom, Sto: undefined")
  }
  


  k<-0
  
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
 
    if( !exists( "paramSpace"  ) | !addToParamSpace  | nrep>1 ){
      # global, holds parameters in each (or best so far?) iterations. Will be a data.frame
      paramSpace <<- c()
    }
    
    subset<-1:nrow( selectedPeaks )
    # usesubsets takes an integer value. Between reps, will randomly select that many lines from selectedPeaks
    if( !is.null( usesubsets ) ){
      if( !is.wholenumber( usesubsets ) | usesubsets < 0 ){stop("Error: usesubsets must be integer")}
      nsubset <- floor( runif( 1, min=usesubsets , max=nrow(selectedPeaks)+1 ) )
      # this tests whether there's a peak to be forces in selectedPeaks.. if so, force it in subset
      if( ncol( selectedPeaks ) == 3 & any( selectedPeaks[, ncol(selectedPeaks) ] > 0 ) ){
        # force the chosen peak in the subset
        subset<- sort( unique( c( which(  selectedPeaks[,3] > 0 ), 
                                  sample( which( selectedPeaks[,3]==0 ) , nsubset-1 ) ) ) )
      } else {
        subset <- sort( sample( 1:nrow( selectedPeaks ), nsubset ) )
      }
      cat( paste( subset, collapse=",") , "\n")
    }
    
    wd<-rep(1,nrow(selectedPeaks[subset,])  ) 
    
    thisMn<<-99999
    # use GenSA, simulated annealing. No need for replicates, it will go to global.
    if( length(optimFct)==1 & optimFct[1]==2 ){
      #require(GenSA)
      if( is.null(control) ){ 
        stop("Error: must specify control. See ?GenSA.") 
        control=list( maxit=500 ) 
      }
      if( is.null( Sn ) ){
        op <- GenSA( par=NULL ,fn=peakProximity,lower=c( Sfrom , lowerF )  ,   
                   upper=c( Sto, upperF ),  control=control , 
                   selectedPeaks=selectedPeaks[subset,]  , npeaks=nrow( selectedPeaks) , 
                   wd=wd, xonly=xonly , notSeenPenalty=notSeenPenalty, ...  ) 
      } else {
        op <- GenSA( par=NULL ,fn=peakProximity,lower=c( lowerF )  ,   
                     upper=c( upperF ),  control=control , Sn=Sn, 
                     selectedPeaks=selectedPeaks[subset,]  , npeaks=nrow( selectedPeaks) , 
                     wd=wd, xonly=xonly , notSeenPenalty=notSeenPenalty, ...  ) 
      }
        
        
      outputlist[[currentrep]] <- op
      outputlist[[currentrep]]$start<-NULL
      outputlist[[currentrep]]$subset<- paste( subset, collapse=",")
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      outputlist[[currentrep]]$paramSpace<-paramSpace
    }
    
    ###
    
    if( length(optimFct)==1 & optimFct[1]==1 ){ # use optim
      if( is.null(control) ){ 
        stop("Error: must specify control. See ?optim.") 
      }
        # starting values  
         start<-startOptim( Sfrom, Sto, nsubcl , forced= ncol( selectedPeaks ) == 3, Sn=Sn )
         
        op<- optim( par=start , fn=peakProximity, selectedPeaks=selectedPeaks[subset,]  ,
                  verbose= T , control=control, npeaks=nrow( selectedPeaks), Sn=Sn,  
                  notSeenPenalty=notSeenPenalty , wd=wd, xonly=xonly,... ) 
     
        
      outputlist[[currentrep]]<-op
      outputlist[[currentrep]]$start<-start
      outputlist[[currentrep]]$subset<- paste( subset, collapse=",")
      outputlist[[currentrep]]$maxc<-maxc
      outputlist[[currentrep]]$maxsubcldiff<-maxsubcldiff
      
    }
    
  } # for( currentrep in 1:nrep)

  # grid search, grid consist of starting values
  # if optimFct > 2 then take that many grid points for each parameters
  # does not need nrep to be defined, the number of rep is determined by number of starting points
  
if( (length(optimFct)==1 & optimFct[1]>2 ) | length(optimFct)>1 ){
  if( length(optimFct)>1 ){
    if( length(optimFct)==1 ){ optimFct<-rep( optimFct[1], length(upperF)+1) }  
    grid<-list( seq( Sfrom, Sto, len=2*optimFct[1]+1 ) )
    grid[[1]]<-grid[[1]][  seq(2,length(grid[[1]]),2)  ]
    
    for( co in 1:length( lowerF)  ){
      grid[[co+1]]<-seq( lowerF[co], upperF[co], len=2*optimFct[co+1]+1 )
      grid[[co+1]]<-grid[[co+1]][ seq(2,length(grid[[co+1]]),2)  ]
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
      op<- optim( par=start , fn=peakProximity, selectedPeaks=selectedPeaks[subset,],  
                  lower=c( Sfrom , lowerF )  ,   upper=c( Sto, upperF ),
                  method=method,
                  verbose= T , control=control, npeaks=nrow( selectedPeaks) , 
                  wd=wd, xonly=xonly , notSeenPenalty=notSeenPenalty , ... ) 
      } else {
        op<- optim( par=start , fn=peakProximity, selectedPeaks=selectedPeaks[subset,],  
                    lower=c( lowerF )  ,   upper=c( upperF ),
                    method=method,
                    verbose= T , control=control, npeaks=nrow( selectedPeaks) , Sn=Sn,
                    wd=wd, xonly=xonly , notSeenPenalty=notSeenPenalty , ... ) 
      }
        
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

## alias 
## surveyParamSpace<-coverParamSpace

