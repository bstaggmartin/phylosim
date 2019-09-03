gen.vec.field<-function(no.dim=2,dim=rep(10,no.dim)){
  stopifnot(length(dim)==no.dim)
  apply(array(rnorm(no.dim*prod(dim)),
              dim=c(no.dim,dim)),
        2:(no.dim+1),function(vec) {vec/sqrt(sum(vec^2))})
}

#add some arg for converting between trait value and grid value?

#interesting emergent phenomenon --> might be more efficient to calc 2D noise slices and stitch them together later,
#rather than do it all in one function; it is lightning-fast for iterating through points in only 2 or fewer dimensions,
#presumably due to precalculation of grid cell coordinates?

#becomes way too slow for >6D noise...

get.noise.value<-function(points,vec.field,interp.fun="cubic"){
  #determine dimensionality of vector field
  if(is.vector(vec.field)==T){
    no.dim<-1
  }else{
    no.dim<-dim(vec.field)[1]
    #if only one interpolation function is specified, extrapolate interpolation function to all dimensions of vector
    #field (converting from vector to character list in process)
    if(length(interp.fun)==1){
      interp.fun<-rep(list(interp.fun),no.dim)
    }
  }
  #convert potential vector of interpolation functions to list
  interp.fun<-as.list(interp.fun)
  #replace any named functions in list of interpolation functions with their actual functions
  interp.fun[interp.fun=="linear"]<-list(function(x) x)
  interp.fun[interp.fun=="cubic"]<-list(function(x) -2*x^3+3*x^2)
  interp.fun[interp.fun=="quintic"]<-list(function(x) 6*x^5-15*x^4+10*x^3)
  #check that all elements of interpolation functions list are functions
  stopifnot(all(unlist(lapply(interp.fun,is.function))))
  if(is.null(nrow(points))){
    no.points<-1
  }else{
    no.points<-nrow(points)
    points<-as.matrix(points)
    colnames(points)<-NULL
  }
  base.grid.point<-floor(points)
  #make sure point and specified interpolation functions are of the same dimensionality as vector field
  stopifnot(no.dim==length(points)/no.points&no.dim==length(interp.fun))
  #find the lowest adjacent grid point
  base.grid.point<-floor(points)
  #make sure lowest adjacent grid point is no less than 1 and not greater than or equal the total number of grid
  #points in any dimension (i.e., specified points lie within grid defined by vector field)
  stopifnot(base.grid.point>=1&base.grid.point<dim(vec.field)[2:(no.dim+1)])
  if(no.points==1){
    #create a matrix of all adjacent grid points by forming all combinations between base grid coordinates and base grid
    #coordinates + 1
    adj.grid.points<-as.matrix(expand.grid(as.list(split(c(base.grid.point,base.grid.point+1),rep(1:no.dim,2)))))
    #create a matrix of grid point vectors
    if(no.dim==1){
      grid.point.vecs<-matrix(vec.field[adj.grid.points],ncol=no.dim,byrow=T)
    }else{
      #create a matrix of indices indicating each element of all adjacent grid point vectors
      indices<-cbind(1:no.dim,adj.grid.points[rep(1:nrow(adj.grid.points),each=no.dim),])
      grid.point.vecs<-matrix(vec.field[indices],ncol=no.dim,byrow=T)
    }
    #create a matrix of vectors from grid points to p
    vecs.to.p<-matrix(rep(points,nrow(adj.grid.points)),ncol=no.dim,byrow=T)-adj.grid.points
    #create a vector of dot rpoducts of grid point vectors and vectors to p
    quads<-2^(no.dim:0)
    dot.products<-rep(0,sum(quads))
    dot.products[1:quads[1]]<-rowSums(grid.point.vecs*vecs.to.p)
    #find interpolated relative points
    #see if there's a way to do this without a for loop
    rel.p<-rep(0,no.dim)
    for(i in 1:no.dim){
      rel.p[i]<-interp.fun[[i]]((points-base.grid.point)[i])
    }
    for(i in 1:no.dim){
      for(j in 1:quads[i+1]){
        dot.products[sum(quads[1:i])+j]<-(1-rel.p[i])*dot.products[sum(c(0,quads)[1:i])+(j*2)-1]+rel.p[i]*dot.products[sum(c(0,quads)[1:i])+j*2]
      }
    }
    dot.products[length(dot.products)]
  }else{
    ###DO THIS FOR ALL UNIQUE FLOOR BASE.GRID.POINT VECTORS
    unique.base.grid.point<-unique(base.grid.point)
    unique.base.grid.point<-split(unique.base.grid.point,1:nrow(unique.base.grid.point))
    unique.adj.grid.points<-lapply(unique.base.grid.point,function(x) as.matrix(expand.grid(as.list(split(c(x,x+1),rep(1:no.dim,2))))))
    if(no.dim==1){
      unique.grid.point.vecs<-lapply(unique.adj.grid.points, function(x) matrix(vec.field[x],ncol=no.dim,byrow=T))
    }else{
      unique.indices<-lapply(unique.adj.grid.points, function(x) cbind(1:no.dim,x[rep(1:nrow(x),each=no.dim),]))
      unique.grid.point.vecs<-lapply(unique.indices, function(x) matrix(vec.field[x],ncol=no.dim,byrow=T))
    }
    noise.values<-rep(0,no.points)
    for(p in 1:no.points){
      ###CALL ADJACENT GRID POINTS FOR BASE GRID POINTS OF CURRENT POINT
      list.index<-which(sapply(unique.base.grid.point,function(x) identical(x,base.grid.point[p,])))
      #create a matrix of vectors from grid points to p
      vecs.to.p<-matrix(rep(points[p,],nrow(unique.adj.grid.points[[list.index]])),ncol=no.dim,byrow=T)-unique.adj.grid.points[[list.index]]
      #create a vector of dot rpoducts of grid point vectors and vectors to p
      quads<-2^(no.dim:0)
      dot.products<-rep(0,sum(quads))
      dot.products[1:quads[1]]<-rowSums(unique.grid.point.vecs[[list.index]]*vecs.to.p)
      #find interpolated relative points
      #see if there's a way to do this without a for loop
      rel.p<-rep(0,no.dim)
      for(i in 1:no.dim){
        rel.p[i]<-interp.fun[[i]]((points[p,]-base.grid.point[p,])[i])
      }
      for(i in 1:no.dim){
        for(j in 1:quads[i+1]){
          dot.products[sum(quads[1:i])+j]<-(1-rel.p[i])*dot.products[sum(c(0,quads)[1:i])+(j*2)-1]+rel.p[i]*dot.products[sum(c(0,quads)[1:i])+j*2]
        }
      }
      noise.values[p]<-dot.products[length(dot.products)]
    }
    noise.values
  }
}