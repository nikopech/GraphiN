#visualise

visualiseGenes<-function(data)
  {udata=unique(data)
  label=match(data,udata)
  list("name"=udata,"label"=label)}  