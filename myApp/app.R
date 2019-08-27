if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()



library(dplyr)
library(shiny) 
library(stringdist)
library(igraph)
library(visNetwork)
library(CINNA)
library(DT)
library(ggraph)
library(graphlayouts)
library(cluster)

source('C:/Users/Paul Kallin/Desktop/myApp/stringdistances.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/myApp/visualisation.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/myApp/filtered.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/myApp/mstClustering.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/myApp/multiple_clustering.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/myApp/visual extended.R', echo=TRUE)
source('C:/Users/Paul Kallin/Desktop/myApp/clusterMetrics.R', echo=TRUE)
# Run the application 
shinyApp(ui = ui, server = server)

