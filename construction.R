graphConstruction<-function(data,sequence=NULL)
  {
  
  data=data[, colSums(is.na(data)) != nrow(data)]
  
  seq=data[,sequence]
  seq<- as.character(seq)
  
  sim=diag(0,length(seq) , length(seq))
  
  for (i in 1:(length(seq)-1)){
    for (j in (i+1):(length(seq))){
      sim[i,j]=stringsim(seq[i],seq[j], method = "osa", useBytes = FALSE)
      sim[j,i]=sim[i,j]
    }}
  
  ig <- graph.adjacency(sim, mode="undirected", weighted=TRUE)
  vertex.attributes(ig)<-as.list(data)
  
  
  plot(ig)
  return (ig)
}
