# stringdistances<-function(seq,algo){
#   
#   
#   sim=diag(0,length(seq) , length(seq))
#   
#   for (i in 1:(length(seq)-1)){
#     for (j in (i+1):(length(seq))){
#       sim[i,j]=1-stringsim(seq[i],seq[j], method = algo, useBytes = FALSE)
#       if (sim[i,j]==0) {sim[i,j]=0.00001}
#       sim[j,i]=sim[i,j]
#     }}
#   
#   return (sim)}


stringdistances<-function(seq,algo){


  sim=diag(0,length(seq) , length(seq))

  for (i in 1:(length(seq)-1)){
    for (j in (i+1):(length(seq))){
      sim[i,j]=1-stringsim(seq[i],seq[j], method = algo, useBytes = FALSE)
      if (sim[i,j]==0) {sim[i,j]=0.00001}
      sim[j,i]=sim[i,j]
    }}

  
  return (as.dist(sim))}