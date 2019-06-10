#wilcoxon


wilcoxon_sum_rank<-function(x,labels,alter)
  {
  
    l=sort(unique(labels))
    pval=diag(0,length(l) , length(l))
    
    for(i in 1:length(l)){
      for (j in 1:length(l))
        {
        x1=x[labels==l[i]]
        x2=x[labels==l[j]]
        pval[i,j]=wilcox.test(x=x1,y=x2,paired=FALSE,alternative =alter )$p.value
      }}
    rownames(pval)=l
    colnames(pval)=l
    return(pval)
}

