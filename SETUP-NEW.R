install.packages(c("stringdist","CINNA","igraph","factoextra","FactoMineR"))
install.packages("graphlayouts")
install.packages("rlang")
install.packages("ggplot2")
install.pavkages("MLmetrics")
install.packages("netrankr")
install.packages("igaph")
install.packages("ggraph")
install.packages("igraphdata")
library("rlang")
library("ggplot2")
library("stringdist")
library("CINNA")
library("igraph")
library("FactoMineR")
library ("factoextra")
library("MLmetrics")
library("netrankr")
library("ggraph")   
library("graphlayouts")
library("igraphdata")


data=read.csv("filterInTableTAG2.txt",sep="",header = TRUE) #data reading

ig=graphConstruction(data,sequence="IMGT.gapped.AA.sequences.V.D.J.REGION") 

#some random graph functions

ig2=delete.edges(ig,which(E(ig)$weight>0.50))
                 
attr=vertex.attributes(ig2)               

plot(ig2)

plot(ig2,vertex.color=visualiseGenes(V(ig)$Summary.V.GENE.and.allele))

#filtering

udata=includeInGraph(data,c("Summary.J.GENE.and.allele","Homsap TRBJ1-2*01 F"))
udata2=excludeFromGraph(udata,c("Summary.D.GENE.and.allele","Homsap TRBD1*01 F"))

ig1=graphConstruction(udata) 
ig2=graphConstruction(udata2)
ig2=delete.edges(ig2,which(E(ig2)$weight>0.6))


#PCA analysis of centrality



#find central nodes

#clsutering
model1=communities_general(ig,threshold=0.6,weights=NULL,algorithm="louvain")

model2=mstClustering(ig)

#Visualise  clusters
a=visualiseGenes(data=data$Summary.D.GENE.and.allele)

#Confusion matrix
mat=confusion(model1$membership,a$label,0.5)





