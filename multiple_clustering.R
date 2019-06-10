#cluster general


communities_general<-function(ig,threshold=1,weights=NULL,algorithm="louvain"){
  
ig=delete_edges(ig,which(E(ig)$weight>threshold))

if (algorithm=="edge_betweenness")
      {if (is.null(weights))
      {model=cluster_edge_betweenness(ig,weights=NULL)}
      else
      {model=cluster_edge_betweenness(ig)}} 
else  
      {if (is.null(weights))
        E(ig)$weight[E(ig)$weight!=0]=1-E(ig)$weight[E(ig)$weight!=0]
         
      
      if (algorithm=="louvain")
        {model=cluster_louvain(ig,weights=weights)}
      else if (algorithm=="fast_greedy")
        {model=cluster_fast_greedy(ig,weights=weights)}
      else if (algorithm=="label_propagation")
        {model=cluster_label_prop(ig,weights=weights,initial= sample(1:30,length(V(ig)),replace=TRUE))}
      else if (algorithm=="walktrap")
        {model=cluster_walktrap(ig,steps=10)}
      else if (algorithm=="leading_eigenvalue")
        {model=cluster_leading_eigen(ig,weights=weights)}
  
      if (is.null(weights))
        E(ig)$weight=1-E(ig)$weight}
 
  plot(ig,vertex.color=(model$membership))
  return(model)
  }
