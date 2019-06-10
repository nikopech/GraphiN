#clustering based on mst


mstClustering<-function(ig,threshold=1/30*dim(ig[])[1]){
  
  mstig=mst(ig,algorithm='prim')
  plot(mstig)
  dis=distances(mstig)
  degree=centr_degree(mstig)
  centroids=which(degree$res>threshold)
  
  mins=apply(dis[,centroids],1,which.min)
  clusters=match(centroids[mins],centroids)
  
  
  plot(mstig,vertex.color=clusters)
  list("clusters"=clusters,"centroids"=centroids)
  }