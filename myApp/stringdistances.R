stringdistances<-function(seq){
  
  
  sim=diag(0,length(seq) , length(seq))
  
  for (i in 1:(length(seq)-1)){
    for (j in (i+1):(length(seq))){
      sim[i,j]=stringsim(seq[i],seq[j], method = "osa", useBytes = FALSE)
      sim[j,i]=sim[i,j]
    }}
  
  return (sim)
  }