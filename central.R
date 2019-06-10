#find central nodes


ig2=delete.edges(ig,which(E(ig)$weight>0.57))
prop=proper_centralities(ig2)

central=calculate_centralities(ig2,include=prop[c(-2,-27,-5,-6,-43,-38:-39,-31:-33)])

central=central[lengths(central)!=0]
x=as.data.frame(central,col.names = names(central))



#ranks=as.data.frame(apply(x,2,order))
central_ranked=as.data.frame(apply(x,2,sort))
colnames(central_ranked)=colnames(x)


#centrality indices needing no reverse order
b=c("Average.Distance","clustering.coefficient","Local.Bridging.Centrality","Wiener.Index.Centrality")

#revert the order of indices 
central_ranked[,-match(b,colnames(central_ranked))]=central_ranked[nrow(central_ranked):1,-match(b,colnames(central_ranked))]




#positions:matrix whose  i-th row contains the position of i-th node in the ranking of every index  
positions=c()
for (i in 1:length(x))
  {positions_first=match(x[,i],central_ranked[,i])
   positions_last=match(x[,i],central_ranked[nrow(x):1,i])
   positions=cbind(positions,(positions_first-positions_last+nrow(x)+1)/2)
   }

colnames(positions)=colnames(x)

#positions.summary
sum=t(apply(positions,1,summary))



#PCA analysis of centrality
model= PCA(x, graph = FALSE)
fviz_screeplot(model, ncp=10)

head(model$var$cos2)

fviz_pca_contrib(model, choice = "var",axes=c(1:5))

fviz_pca_var(model, col.var="cos2") +scale_color_gradient2(low="white", mid="blue", high="red", midpoint=0.5) + theme_minimal()

fviz_pca_var(model, col.var="contrib") +scale_color_gradient2(low="white", mid="blue", high="red", midpoint=1.5) + theme_minimal()

fviz_pca_contrib(model, choice = "ind", axes = c(1:1))

fviz_pca_ind(model, col.ind="contrib")+scale_color_gradient2(low="white", mid="blue", high="red", midpoint=1) + theme_minimal()


#clustering-centralities combo

model1=communities_general(ig,threshold=0.6,weights=NULL,algorithm="louvain")
compared_clusters=wilcoxon_sum_rank(x=x[,"Closeness.Centrality..Freeman."],labels=model1$membership,"greater")

