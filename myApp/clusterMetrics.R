#clustering metrics 


conductance<-function(graph,membership)
{
  
  edges=as_data_frame(graph)
  edges$from=membership[edges$from]
  edges$to=membership[edges$to]
  con=c()
  total=sum(edges$weight)
  for (i in sort(unique(membership))){     
            inter=sum(edges$weight[xor(edges$to==i,edges$from==i)])
            clust=sum(edges$weight[edges$to==i | edges$from==i])
            con=c(con,inter/min(total-clust,clust))
  }
  list (conductance=1-mean(con),conductances=con)
}


coverage<-function(graph,membership)
{
  edges=as_data_frame(graph)
  edges$from=membership[edges$from]
  edges$to=membership[edges$to]
  return(sum(edges$weight[edges$from==edges$to])/sum(edges$weight))
}
  