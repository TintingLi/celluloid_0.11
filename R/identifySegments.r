

identifySegments<-function( x,y,segments ){
 id<-identify( x,y, n=1 )
 while( length(id)>0 ){
   print( segments[id,] )
    id<-identify( x,y, n=1 )
 }
}


  


