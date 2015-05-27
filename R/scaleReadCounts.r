

scaleReadCounts <-function( tumourrangedata , epp  ){

 # this function scales tc (the irange object)

 # epp contains the expected peaks, as well as "scaled read counts(say) " 
 # as x and allelic ratio as ar.
 # m1, p1, m2, p2, etc  are counts of maternal (say) and paternal alleles in tumour subclones 1 and 2 

 # for two subclone system:
 # names(epp):
 # [1] "m0" "p0" "m1" "p1" "m2" "p2" "x"  "ar"


 # I want to transform the HMM corrected copy count  ($copy object) 

 mat<-epp[,seq( 3, ncol(epp)-2,2)] 
 pat<-epp[,seq( 4, ncol(epp)-2,2)] 

# DEBUGGED in v7, added data.frame in case vector
 matpat<-data.frame(mat+pat)

 sel <- apply( matpat==0, 1, all )
 zero<- epp[sel, "x" ][1]

 sel <-  apply( matpat==1 , 1, all )
 one<-  epp[sel,"x" ][1]
 
 sel<- apply( matpat==2 , 1, all )
 two<- epp[sel,"x" ][1]

 # this preserves the ratios 


 cntransform<-function( x ){ (x-zero)/(two-zero)   } 
 
# NOTE: log is not appropriate scale becuase of normal contamination.
# it's not because read count is reduced by half that the copy number is reduced by half as well 

 tumourrangedata$icopy <- 2*cntransform( 2^tumourrangedata$copy ) 

# sel1<- tumourrangedata$icopy > 0 & !is.na(tumourrangedata$icopy) 
# sel2<- tumourrangedata$icopy <= 0 & !is.na(tumourrangedata$icopy)
# tumourrangedata$icopy[sel1]<- log( tumourrangedata$icopy[sel1], 2 )
# tumourrangedata$icopy[sel2]<- NA 

 return(tumourrangedata) 

}


scaleSegments <-function( seg , epp  ){

 seg$imean<-seg$mean

 mat<-epp[,seq( 3, ncol(epp)-2,2)] 
 pat<-epp[,seq( 4, ncol(epp)-2,2)] 

 matpat<-data.frame(mat+pat)

 sel <- apply( matpat==0, 1, all )
 zero<- epp[sel, "x" ][1]

 sel <-  apply( matpat==1 , 1, all )
 one<-  epp[sel,"x" ][1]
 
 sel<- apply( matpat==2 , 1, all )
 two<- epp[sel,"x" ][1]

 cntransform<-function( x ){ (x-zero)/(two-zero)   } 
 
# NOTE: log is not appropriate scale becuase of normal contamination.
# it's not because read count is reduced by half that the copy number is reduced by half as well 

 seg$imean  <- 2*cntransform( seg$mean ) 

 return(seg) 

}





