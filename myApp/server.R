#tempdata=data[data[,l[1]] %in% l[2:length(l)]]

# use the below options code if you wish to increase the file input limit, in this example file input limit is increased from 5MB to 9MB
# options(shiny.maxRequestSize = 9*1024^2)

shinyServer(function(input,output,session){
  
  values <- reactiveValues(flagEX=0,filterdf=data.frame(Columns=character(),Keys=character(),"I/E"=character()),ig=NULL,indexes=c(),forest=c(),sim=NULL)

  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
    fulldata <- eventReactive(input$file,{
        file1 <- input$file
        if(is.null(file1)){return()} 
        data=read.csv(file=file1$datapath, sep=input$sep, header = input$header)
        data=data[sample(max(dim(data)),300, replace = FALSE),]
        #sim=stringdistances(seq=as.character(data[,input$seqSelect]),algo=input$simSelect)
        values$indexes=1:nrow(data)
        data
       })
    
    
    observeEvent(input$simButton,{
            
            print(system.time({sim=stringdistances(seq=as.character(fulldata()[,input$seqSelect]),algo=input$simSelect)}))
            shinyalert("Distances Calculated", type = "success")
            values$sim=sim
      
      })
  
    
    observeEvent(input$graphButton, {
         if (is.null(values$sim )) {return()}
         templist=list(fulldata()[values$indexes,],as.matrix(values$sim)[values$indexes,values$indexes])
          
         if (input$clusterId)
          { 
               uni=unique(templist[[1]][,c(input$seqSelect,"cluster_id")])
               ind=1
               for (i in 2:nrow(uni))
                  { floorindex=ind[i-1]+1
                    ind=c(ind,ind[i-1]+min(which(templist[[1]][floorindex:nrow(templist[[1]]),input$seqSelect]==uni[i,1] & templist[[1]][floorindex:nrow(templist[[1]]),"cluster_id"]==uni[i,2] )))
                   }
               templist[[1]]=templist[[1]][ind,]
               templist[[2]]=templist[[2]][ind,ind]
          }
          
         
          
          ig<- graph.adjacency(templist[[2]], mode="undirected", weighted=TRUE)
          ig=delete.edges(ig,which(E(ig)$weight>input$slider))
          vertex.attributes(ig)<-as.list(templist[[1]])
          
          V(ig)$value=V(ig)$freq_cluster_id^(1/3.5)/10
          
          V(ig)$color.background=colors[visualiseGenes(data=V(ig)$dataName)$label]
          V(ig)$color.border=V(ig)$color.background
        
          
          m=(1-E(ig)$weight-min(1-E(ig)$weight))/(max(1-E(ig)$weight)-min(1-E(ig)$weight))
          
          
          E(ig)$width <- m*2.5
          E(ig)$title<-E(ig)$weight
          
          V(ig)$id=1:gorder(ig)
          V(ig)$label<-rownames(templist[[1]])
          forest=components(ig,mode="weak")
          
          x=forest[[1]]
          y=table(x)
          z=y[match(x,names(y))]
          x[z<5]=NA
          x=match(x,sort(unique(x),na.last = TRUE))
          values$forest=x
          values$ig=ig
          choices=as.list(sort(values$forest))
          names(choices)=sort(values$forest)
          updateSelectInput(session,"componentSelect",choices=choices,selected=1)
          clusterValues$member=NULL
          clusterValues$membership=vector(mode="list",length=7)
          clusterValues$comm1=NULL
          clusterValues$comm2=NULL
          
          names(clusterValues$membership)=c("louvain","fast_greedy","label_propagation","leading_eigenvalue","walktrap","edge_betweenness","hierarchical")
          updateSelectInput(session,"clusterSelect",selected=" ")
          mstValues$edges=NULL
          
          })
    
    
    observeEvent(input$componentButton,{
          if (is.null(values$ig)) {return()}
          values$ig=induced_subgraph(values$ig,which(values$forest==input$componentSelect),impl="auto")
          V(values$ig)$id=1:gorder(values$ig)
          values$forest=c()
          values$forest[1:gorder(values$ig)]=1
          updateSelectInput(session,"componentSelect",choices=list("1"=1),selected=1)
          clusterValues$member=NULL
          clusterValues$membership=vector(mode="list",length=7)
          clusterValues$comm1=NULL
          clusterValues$comm2=NULL
          names(clusterValues$membership)=c("louvain","fast_greedy","label_propagation","leading_eigenvalue","walktrap","edge_betweenness","hierarchical")
          updateSelectInput(session,"clusterSelect",selected=" ")
          mstValues$edges=NULL
       })
  
        
    output$dataset<-DT::renderDataTable({
          temp=fulldata()
          if (is.null(temp)){return()}
          if (input$idInput!=0)
             {temp=temp[input$idInput==rownames(temp),]}
          temp
          
          })
    
    
    
    output$network <- renderVisNetwork({
          xx=input$visGraph
          igtemp=isolate(values$ig)
          if (is.null(igtemp)) {return()}
          
          coords=layout_with_stress(igtemp)
            
          data_vertices=as_data_frame(igtemp,what=c("vertices"))[,c("value","label","color.border","color.background","id")]
          data_vertices=cbind(data_vertices,title=isolate(values$forest))
          data_edges=as_data_frame(igtemp,what=c("edges"))
          graph<-visNetwork(data_vertices, data_edges) %>%
            visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
            visOptions (nodesIdSelection = list("useLabels"=TRUE),selectedBy = list("variable"="title"),highlightNearest = TRUE)   %>%
            #visPhysics(solver="repulsion") %>%
            visInteraction(multiselect = TRUE)%>%
            visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) %>%
            #visIgraphLayout(layout = "layout_with_fr",smooth=FALSE,physics = FALSE) %>%
            visLegend(addNodes = data.frame()) 
          
          graph
         
          
          
         
        })
          
          
          
    output$selectbox1<-renderUI({
          tempdata=fulldata()
         
          choice=list()
          for (l in colnames(tempdata))
             choice=c(choice,l)
          names(choice)=colnames(tempdata)
          selectInput("select1", label = "Select a column", choices = choice)
        })
       
       
       
       
         
    output$textbox1<-renderUI({
          tempdata=fulldata()
          selected=input$select1
          if (is.null(selected)) {return()}
          if (is.numeric(tempdata[,selected]))
             {sliderInput("slider1", label = "", min = min(tempdata[,selected]), max = max(tempdata[,selected]), value = c(min(tempdata[,selected]),max(tempdata[,selected]) ))}
          else
              {textInput("text1","")}
        
        }) 
     
 
      
      

      
      xxchange <- reactive({
        paste(input$excludeButton ,input$includeButton)
      })
     
      
      
      observeEvent(xxchange(),
                      {
                     if(input$excludeButton[[1]]==0 & input$includeButton[[1]]==0) {return()}
                     tempdata=fulldata()  
                                
                     if (isolate(values$flagEX)!=input$excludeButton[[1]]){
                                    ie="E"
                                    values$flagEX=input$excludeButton[[1]]}
                     else {ie="I"}
                                
                     selected=input$select1
                     if (is.numeric(tempdata[,selected])){ 
                                  x=data.frame("Columns"=selected,"Keys"=paste0(input$slider1[1],",",input$slider1[2]),"I/E"=ie)}
                     else{
                                  x=data.frame("Columns"=selected,"Keys"=input$text1,"I/E"=ie)}
                                
                   
                   values$filterdf<-rbind(values$filterdf,x)
                                
                              })
      
      
      observeEvent(input$ResetButton,{
                     values$filterdf<-data.frame(Columns=character(),Keys=character(),"I/E"=character())
                    })
      
      
      
      observeEvent(input$FilterButton,{
                  if (nrow(values$filterdf)==0) {return()}
                  tempdata=isolate(fulldata())
                  filterg=values$filterdf
                  exc=filterg[filterg[,3]=="E",1:2]
                  inc=filterg[filterg[,3]=="I",1:2]
                  filterinc=list()
                  filterexc=list()
                  
                  if ((nrow(exc))!=0){
                    for (i in 1:nrow(exc)){
                      filterexc[[i]]=c(as.character(exc[i,1]),unlist(strsplit(as.character(exc[i,2]),",")))
                      }}
                    
                  if ((nrow(inc))!=0){
                    for (i in 1:nrow(inc)){
                      filterinc[[i]]=c(as.character(inc[i,1]),unlist(strsplit(as.character(inc[i,2]),",")))
                    }}  
                      
                      
                  pointers1=includeInGraph(tempdata,filterinc)
                  if (length(pointers1)!=0) 
                    {pointers2=excludeFromGraph(tempdata[pointers1,],filterexc)
                     if (length(pointers2)!=0){
                       values$indexes=pointers1[pointers2]}
                     else {values$indexes=NULL}}
                  else
                    {values$indexes=NULL}
                  
                  if (input$ReverseButton==TRUE){values$indexes=setdiff(1:nrow(fulldata()),values$indexes)}
                 },ignoreInit = TRUE)                 
      
      
      
      output$Indexes<-renderText({
            if (is.null(values$indexes)) {return()}
            values$indexes
      })
      
      
      
      output$filtertable<-renderTable({
                          
            if (is.null(values$filterdf)){return()}
            values$filterdf 
       })
      
      
      output$downloadGraph <- downloadHandler(
          filename = function() {
              paste("graph.zip")
            },
          content = function(file) {
              fileName<- paste("edges.csv")
              write.csv(as_data_frame(values$ig)[,c("from","to","weight")], fileName, row.names = FALSE)
              files=fileName
              fileName<-paste("vertices.csv")
              write.csv(V(values$ig)$label, fileName, row.names = FALSE)
              files<-c(fileName,files)
              zip(file,files)
              
      })
      
      
      
      
      
      
      mstValues<-reactiveValues(edges=NULL,centroids=NULL,clusters=NULL,flag2=NULL,legenddf=NULL,legenddf2=NULL,shape=NULL)
      
      
      # observeEvent(input$mstButton,{
      # 
      #         if (is.null(values$ig)) {return()}
      #         mstig=mst(values$ig,algorithm='prim')
      #         a=mstClustering(values$ig)
      # 
      # 
      # 
      # 
      #         V(mstig)$value=V(values$ig)$freq_cluster_id^(1/3.5)*1
      # 
      #         m=(1- E(mstig)$weight-min(1-E(mstig)$weight))/(max(1-E(mstig)$weight)-min(1-E(mstig)$weight))
      #         E(mstig)$width <- m*2.5
      #         n=c()
      #         n[gorder(mstig)+1]=1
      #         n[a]="triangle"
      #         V(mstig)$shape=n[1:gorder(mstig)]
      #         V(mstig)$id=1:gorder(mstig)
      #         mstValues$ig=list(vertices=as_data_frame(mstig,what="vertices")[,c("id","label","shape","value")],edges=as_data_frame(mstig))
      #         mstValues$centroids=V(mstig)$label[a] #a$centroids]
      #         mstValues$flag2=NULL
      # 
      # },ignoreInit = TRUE)

      
      observeEvent(input$mstButton,{

            if (is.null(values$ig)) {return()}
            igtemp=values$ig
            if (input$mstAlgoSelect=="Kruskal")
                  edges=MST(igtemp,input$mstAlgoSelect)
            else
                  edges=as_data_frame(mst(igtemp,algorithm = "prim"))
            
            #a=mstClustering(values$ig)
            ig=graph_from_data_frame(edges,directed=FALSE,vertices=as_data_frame(igtemp,what="vertices")[,c("id","label","value")])
            
            degree=centr_degree(ig)
            centroids=which(degree$res>1/15*gorder(ig))
            
            
            n=c()
            n[gorder(ig)+1]=1
            n[centroids]="triangle"
            mstValues$shape=n[1:gorder(ig)]
            



            #V(mstig)$value=V(values$ig)$freq_cluster_id^(1/3.5)*1

            m=(1-edges$weight-min(1-edges$weight))/(max(1-edges$weight)-min(1-edges$weight))
            edges$width <- m*2.5
           
            #V(mstig)$id=1:gorder(mstig)
            mstValues$edges=edges
            mstValues$centroids=V(igtemp)$label[centroids] #a$centroids]
            mstValues$flag2=NULL

      },ignoreInit = TRUE)
      
      
      output$mstnetwork<-renderVisNetwork({
            xx=input$visMST
            isolate({
            igtemp=values$ig
            edges=mstValues$edges
            if (is.null(edges) | is.null(igtemp)) {return()}
            
            ig=graph_from_data_frame(edges,directed=FALSE,vertices=as_data_frame(igtemp,what="vertices")[,c("id","label","value")])
            
            degree=centr_degree(ig)
            centroids=which(degree$res>1/15*gorder(ig))
            coords=layout_with_stress(ig)
            #coords=layout_with_kk(igtemp,coords = coords)
            
            V(ig)$id=V(ig)$name
            delete_vertex_attr(ig,"name")
            
            
                  if ((input$clusterMST) & !is.null(clusterValues$member) ){   
                        cbg=colors[clusterValues$member]
                        mstValues$legenddf=data.frame("label"=unique(clusterValues$member),"color"=(colors[unique(clusterValues$member)]))}
                  
                  else{  
                        if (input$colormst=="Default"){
                            cbg="black"
                            mstValues$legenddf=NULL}
                        else
                            {temp=visualiseGenes(get.vertex.attribute(igtemp,input$colormst))
                             cbg=colors[temp$label]
                             mstValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}}#unique(temp$label)
                  
                  if (input$bordermst=="Default")
                        {cbor="black" #colors[a$clusters]
                        mstValues$legenddf2=NULL}
                  else
                        {temp=visualiseGenes(get.vertex.attribute(igtemp,input$bordermst))
                        cbor=colors[temp$label]
                        mstValues$legenddf2=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}
            
           
            
            
            mstValues$flag2=c(1)
            visNetwork(cbind(as_data_frame(ig,what="vertices"),"color.background"=cbg,"color.border"=cbor,"shape"=mstValues$shape),cbind(as_data_frame(ig),"title"=E(ig)$weight)) %>%
                visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
                visNodes(borderWidth = 2) %>%
                visOptions (nodesIdSelection = TRUE,highlightNearest = TRUE)   %>%
                ##visPhysics(solver="forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant=-2000,centralGravity=0.02,springLength=40,springConstant=0.4,damping=1,avoidOverlap=1)) %>%
                visInteraction(multiselect = TRUE)%>%
                visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) 
            })
                
       })
      
      
      output$mstLegend<-renderPlot({
            if (is.null(mstValues$legenddf)){return()}
            mstValues$legenddf=mstValues$legenddf[order(mstValues$legenddf$label),]
            plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
            legend("bottomleft", legend =mstValues$legenddf$label, pch=16, pt.cex=3, cex=1.5, bty='n',
                   col = as.character(mstValues$legenddf$color))
            mtext("Background",side=2, at=0.2, cex=2)
        
      })
      
      
      
      output$mstLegend2<-renderPlot({
            if (is.null(mstValues$legenddf2)){return()}
            mstValues$legenddf2=mstValues$legenddf2[order(mstValues$legenddf2$label),]
            plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
            legend("bottomleft", legend =mstValues$legenddf2$label, pch=16, pt.cex=3, cex=1.5, bty='n',
                   col = as.character(mstValues$legenddf2$color))
            mtext("Border",side=2, at=0.2, cex=2)
        
      })
      
      output$Cendroids<-renderText({
            if (is.null(mstValues$centroids)) {return()}
            mstValues$centroids
      })
      
      
      output$Links<-renderText({
            if (is.null(mstValues$keyvertices)) {return()}
            mstValues$keyvertices
      })
     
      
      
      observe({
            tempdata=fulldata()
            if (is.null(tempdata)) {return()}
            choice=list("Default")
    
            for (l in colnames(tempdata)){
              if (!is.numeric(tempdata[,l])){
                choice=c(choice,l)
                
              }}
            names(choice)=c("Default",unlist(choice)[-1])
            updateSelectInput(session,"colormst", label = "Select an attribute  for nodes background coloring", choices = choice)
            updateSelectInput(session,"bordermst", label = "Select an attribute for nodes border coloring", choices = choice)
            updateSelectInput(session,"clusterColor", label = "Select an  attribute for nodes coloring", choices = choice)
            updateSelectInput(session,"centralColor", label = "Select an  attribute for nodes coloring ", choices = choice[-1],selected="dataName")
            updateSelectInput(session,"clusterSelect2",label="Select a clustering algorithm or a nodes attribute",choices=c(choice[-1],list("Louvain"="louvain",
                                                                                                                               "Fast Greedy"="fast_greedy",
                                                                                                                               "Label Propagation"="label_propagation",
                                                                                                                               "Leading Eigenvalue"="leading_eigenvalue",
                                                                                                                               "Walktrap"="walktrap",
                                                                                                                               "Edge Betweeness"="edge_betweenness")),selected = "louvain")
      })
     
      
        
        observe({
          
              x=input$colormst
              y=input$bordermst
              z=input$clusterMST
            
              
              if ( is.null(isolate(mstValues$flag2)) | is.null(isolate({values$ig}))){return()}
              tempdata=as_data_frame(isolate({values$ig}),what="vertices")
             
              
              
              
              if (!is.null(clusterValues$member) & input$clusterMST)
                  {bgcolor=colors[clusterValues$member]
                   mstValues$legenddf=data.frame("label"=unique(clusterValues$member),"color"=colors[unique(clusterValues$member)])}
                 
                
              else{  
                 if (x=="Default"){
                      bgcolor="black"#colors[mstValues$clusters]
                      mstValues$legenddf=NULL}#data.frame("label"=unique(mstValues$clusters),"color"=unique(colors[mstValues$clusters]))}
                 else {
                      temp=visualiseGenes(tempdata[,input$colormst])
                      bgcolor=colors[temp$label]
                      mstValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}}
              
              if (y=="Default")
                {bordcolor="black" #colors[mstValues$clusters]
                 mstValues$legenddf2=NULL}#data.frame("label"=unique(mstValues$clusters),"color"=unique(colors[mstValues$clusters]))}
              else  
                {temp=visualiseGenes(tempdata[,input$bordermst])
                 bordcolor=colors[temp$label]
                 mstValues$legenddf2=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}
                
              visNetworkProxy("mstnetwork") %>%
                  visUpdateNodes(data.frame(id=1:nrow(tempdata),color.background=bgcolor,color.border=bordcolor))
          
        })
        
        
        # output$downloadMST<- downloadHandler(
        #   filename = function() {
        #     paste("mst.zip")
        #   },
        #   content = function(file) {
        #     fileName<- paste("mst-edges.csv")
        #     write.csv(mstValues$ig[["edges"]][,c("from","to","weight")], fileName, row.names = FALSE)
        #     files=fileName
        #     fileName<-paste("mst-vertices.csv")
        #     write.csv(mstValues$ig[["vertices"]][,"label"], fileName, row.names = FALSE)
        #     files<-c(fileName,files)
        #     zip(file,files)
        #     
        #   })

         output$downloadMST <- downloadHandler(
           filename = function() {
             paste("mst.csv")
           },
           content = function(file) {
             write.csv(mstValues$ig[["edges"]][,c("from","to")],file, row.names = FALSE)
          }
         )
        
        
        
       
        
        
         centralValues=reactiveValues(centralities=data.frame(),rankings=data.frame(),summary=data.frame(),legenddf=NULL)
        
       
         observeEvent(input$centralButton,{
                if (is.null(values$ig)){return()}
                igtemp=giant_component_extract(values$ig)[[1]]
                
                prop=proper_centralities(igtemp)
                print("centralities")
                if (input$fastCButton){
                    print(system.time({
                      central=calculate_centralities(igtemp,include=prop[c(3,10,26)])
                      E(igtemp)$weight=1-E(igtemp)$weight
                      central=c(central,calculate_centralities(igtemp,include=prop[15]))
                      }))
                    b=c("Average.Distance")}
                else { 
                    print(system.time({
                      central=calculate_centralities(igtemp,include=prop[c(-2,-5,-6,-8,-15,-19:-20,-21,-24,-27,-31:-33,-38:-39,-43,-50)])
                      E(igtemp)$weight=1-E(igtemp)$weight
                      central=c(central,calculate_centralities(igtemp,include=prop[c(15,19,20,24,50)]))
                      }))
                    b=c("Average.Distance","Local.Bridging.Centrality","Wiener.Index.Centrality")}
               
                central=central[lengths(central)!=0]
                x=as.data.frame(central,col.names = names(central),row.names=V(igtemp)$label)
               
                print("ranks")  
                print(system.time({
                central_ranked=as.data.frame(apply(x,2,sort))
                colnames(central_ranked)=colnames(x)
                
                
                #revert the order of indices 
                central_ranked[,-match(b,colnames(central_ranked))]=central_ranked[nrow(central_ranked):1,-match(b,colnames(central_ranked))]
              
                positions=c()
                for (i in 1:length(x))
                  { positions_first=match(x[,i],central_ranked[,i])
                    positions_last=match(x[,i],central_ranked[nrow(x):1,i])
                    positions=cbind(positions,(positions_first-positions_last+nrow(x)+1)/2)
                }
                }))
                colnames(positions)=colnames(x)
                rownames(positions)=rownames(x)
                
                centralValues$centralities=x
                centralValues$rankings=positions
                centralValues$summary=as.data.frame(t(apply(positions,1,summary)),rownames=rownames(x))
            
            
          
        },ignoreInit = TRUE)
        
        output$Centralities<-DT::renderDataTable({
                if (is.null(centralValues$centralities)) {return()}
                centralValues$centralities
      
          
          })
       
         output$Rankings<-DT::renderDataTable({
                if (is.null(centralValues$rankings)) {return()}
                centralValues$rankings
          
          
         },rownames=FALSE)
      
         
         output$Summary<-DT::renderDataTable({
           
                 if (is.null(centralValues$summary)) {return()}
                 centralValues$summary
           
          },rownames=FALSE)
         
         output$centralnetwork<-renderVisNetwork({
           
                  igtemp=giant_component_extract(isolate(values$ig))[[1]]
                  if (is.null(igtemp) | is.null(centralValues$centralities)) {return()}
                  
                  
                  bc <- centralValues$centralities[,input$centralSelect]
                  if (input$centralSelect=="Average.Distance") {bc=1/bc}
                  
                  coords=layout_with_centrality(igtemp,cent=bc,iter=500,tol=1e-04)
                  # p1 <- ggraph(igtemp,layout = "centrality", cent = bc,iter = 500, tol = 1e-04) +
                  #   draw_circle(use = "cent") +
                  #   annotate_circle(bc,format="",pos="bottom") +
                  #   geom_edge_link(edge_color="black",edge_width=0.3)+
                  #   geom_node_point(aes(fill=as.factor(Faction)),size=2,shape=21)+
                  #   scale_fill_manual(values=c("#8B2323", "#EEAD0E"))+
                  #   theme_graph()+
                  #   theme(legend.position = "none")+
                  #   coord_fixed()+
                  #   labs(title="betweenness centrality")
                  # 
                  # coords=cbind(x=p1$data$x,y=p1$data$y)
                  
                 
                  
                  i=c()
                  i[gorder(igtemp)+1]=""
                  i[coords[,1]==0 & coords[,2]==0]="triangle"
                  
               
                  
                  temp=visualiseGenes(get.vertex.attribute(igtemp,"dataName"))
                  centralValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])
                  updateSelectInput(session,"centralColor",selected="dataName")
                  updateCheckboxInput(session,"clusterCentral",value=FALSE)
                  
                  data_vertices=cbind(as_data_frame(igtemp,what=c("vertices"))[,c("value","label")],shape=i[1:gorder(igtemp)],color=colors[temp$label])
                  data_vertices$id=rownames(data_vertices)
                  data_edges=data.frame(from=1,to=2,weight=1)
                  
                  visNetwork(data_vertices, data_edges) %>%
                    visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
                    visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
                    #visPhysics(solver="repulsion") %>%
                    visInteraction(multiselect = TRUE)%>%
                    visEdges(color="white")%>%
                    visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) 
          })
         
         
         
         
         observe({
           x=input$clusterCentral
           y=input$centralColor
           
           if (is.null(isolate(values$ig))){return()}
           igtemp=giant_component_extract(isolate(values$ig))[[1]]
           
           if (!is.null(clusterValues$member) && input$clusterCentral)
              {tempcolor=colors[clusterValues$member[get.vertex.attribute(igtemp,"id")]]
               centralValues$legenddf=NULL}
           
            else {
              temp=visualiseGenes(get.vertex.attribute(igtemp,input$centralColor))
              tempcolor=colors[temp$label]
              centralValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}
           
           
           
           visNetworkProxy("centralnetwork") %>%
             visUpdateNodes(data.frame(id=1:gorder(igtemp),color=tempcolor))
           
         })
         
         output$centralLegend<-renderPlot({
           if (is.null(centralValues$legenddf)){return()}
           centralValues$legenddf=centralValues$legenddf[order(centralValues$legenddf$label),]
           plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
           legend("bottomleft", legend =centralValues$legenddf$label, pch=16, pt.cex=3, cex=1.5, bty='n',
                  col = as.character(centralValues$legenddf$color))
           mtext("Nodes Color",side=2, at=0.2, cex=2)
           
         })
         
         
         output$downloadCentral <- downloadHandler(
           filename = function() {
             paste("centralities.zip")
           },
           content = function(file) {
             fileName<- paste("centralities.csv")
             write.csv(centralValues$centralities, fileName, row.names = FALSE)
             files=fileName
             fileName<-paste("rankings.csv")
             write.csv(centralValues$rankings, fileName, row.names = FALSE)
             files<-c(fileName,files)
             fileName<-paste("summary.csv")
             write.csv(centralValues$summary, fileName, row.names = FALSE)
             files<-c(fileName,files)
             zip(file,files)
             
           })
         
         
         
         
         
         
         #centralities=data.frame(),rankings=data.frame(),summary=data.frame()
         
         clusterValues=reactiveValues(member=NULL,comm1=NULL,comm2=NULL,legenddf=NULL,membership=NULL)
         
         
         observeEvent(input$clusterSelect,{
           igtemp=values$ig
           if (is.null(igtemp) ){return()}
           if (input$clusterSelect==" ") {clusterValues$member=NULL}
           else 
                   {    
                   if (is.null(clusterValues$membership[[input$clusterSelect]]))
                       {id <<- showNotification(paste("Clustering Calculation..."), duration = NULL)
                       print("cluster")
                       print(system.time({model=communities_general(igtemp,algorithm =input$clusterSelect)}))
                       clusterValues$membership[[input$clusterSelect]]=model$membership
                       if (!is.null(id)){
                            removeNotification(id)
                            id <<- NULL}
                       }
                  clusterValues$member=clusterValues$membership[[input$clusterSelect]]
                    }
         })
           
         
           
         output$clusterNetwork<-renderVisNetwork({
               
               xx=input$visCluster
               igtemp=isolate(values$ig)
               membertemp=isolate(clusterValues$member)
               if (is.null(membertemp) | is.null(igtemp)) {return()}
               updateSelectInput(session,"clusterColor",selected="Default")
               
               V(igtemp)$grp=membertemp
               print("layout")
               print(system.time({bb <- layout_as_backbone(igtemp,keep=0.4)}))
               V(igtemp)$color=colors[membertemp]
               
               data_vertices=cbind(as_data_frame(igtemp,what=c("vertices")))[,c("value","label","color")]
               data_vertices$id=rownames(data_vertices)
               data_edges=as_data_frame(igtemp,what=c("edges"))
               
             
               
               visNetwork(data_vertices, data_edges) %>%
                 visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
                 visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
                 #visPhysics(solver="repulsion") %>%
                 visInteraction(multiselect = TRUE)%>%
                 visIgraphLayout(layout = "layout.norm", layoutMatrix = cbind(x=bb$xy[,1],y=bb$xy[,2]),smooth=FALSE,physics = FALSE) 
           })
         
         
         
         observeEvent(paste(input$clusterColor,clusterValues$member),{
             if (is.null(values$ig)){return()}
             igtemp=values$ig
             if (input$clusterColor!="Default"){
                  temp=visualiseGenes(get.vertex.attribute(igtemp,input$clusterColor))
                  tempcolor=colors[temp$label]
                  clusterValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[1:length(unique(temp$name))])}
             else {
                  if (!is.null(clusterValues$member))
                     {tempcolor=colors[clusterValues$member]}
                  else
                     {tempcolor=NULL}
                      clusterValues$legenddf=NULL
             }
             if (is.null(tempcolor)){return()}
             visNetworkProxy("clusterNetwork") %>%
               visUpdateNodes(data.frame(id=1:gorder(igtemp),color=tempcolor))
           
          })
         
         
         output$clusterLegend<-renderPlot({
           if (is.null(clusterValues$legenddf)){return()}
           clusterValues$legenddf=clusterValues$legenddf[order(clusterValues$legenddf$label),]
           plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
           legend("bottomleft", legend =clusterValues$legenddf$label, pch=16, pt.cex=3, cex=1.5, bty='n',
                  col = as.character(clusterValues$legenddf$color))
           mtext("Nodes Color",side=2, at=0.2, cex=2)
           
         })
         
         output$metrics<-renderText({
           if (is.null(clusterValues$member)) {return()}
           igtemp=isolate(values$ig)
           E(igtemp)$weight=1-E(igtemp)$weight
           mod=modularity(igtemp,clusterValues$member)
           con=conductance(igtemp,clusterValues$member)$conductance
           cov=coverage(igtemp,clusterValues$member)
           paste0("Modularity:",mod,"  Conductance:",con,"  Coverage:",cov)
           
           
         })
         
         
         
         
         output$intraHeatmap<-renderPlot({
                  if (is.null(clusterValues$member)) {return()}
                  dis=as.matrix(isolate(values$ig[]))
                  maxt=max(apply(dis,1,max))
                  dis[dis==0]=maxt+0.1
                  diag(dis)=0
                  #dis=distances(isolate(values$ig))
                  temp=order(clusterValues$member)
                  dis=dis[temp,temp]
                  heatmap.2(dis,col=colorpanel(256,"red","orange","yellow"),scale="none",trace = "none", density.info = "none",dendrogram="none",Rowv = NA, Colv = NA)
                  #heatmap(dis, Rowv = NA, Colv = NA, col = heat.colors(256), revC = TRUE)
           })
         
         
         
         output$silhouette<-renderPlot({
                 if (is.null(clusterValues$member)) {return()}
                 dis=as.matrix(isolate(values$ig[]))
                 maxt=max(apply(dis,1,max))
                 dis[dis==0]=maxt+0.1
                 diag(dis)=0
                 #dis=distances(isolate(values$ig))
                 plot(silhouette(dist=dis,x=clusterValues$member))
         })
         
         
         
         
         observe({
                igtemp=values$ig
                if (is.null(igtemp)){return()}
               
                 clusterValues$comm1=clusterValues$membership[[input$clusterSelect1]]#communities_general(igtemp,algorithm =input$clusterSelect1)$membership
    
                if (input$clusterSelect2 %in% c("louvain","fast_greedy","label_propagation","leading_eigenvalue","walktrap","edge_betweenness")){
                     
                     clusterValues$comm2=clusterValues$membership[[input$clusterSelect2]]}#communities_general(igtemp,algorithm =input$clusterSelect2)$membership}
                else {
                     clusterValues$comm2=get.vertex.attribute(igtemp,input$clusterSelect2)}
                
                
                #validate(
                #  need(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)), 'Check at least one letter!')
                #)
            })
         
         
         output$confusionMatrix1<-DT::renderDataTable({
                #if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
           
           validate(
             need(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)), "One of the clusterings hasn't been calculated!")
                 )
           
           mat=as.data.frame(t(confusion(clusterValues$comm1,clusterValues$comm2)$mat1))
           mat[,c(1,2)]<-mat[,c(2,1)]
           isolate({colnames(mat)<-c(input$clusterSelect1,input$clusterSelect2,"Freq")})
           mat
           
                })
         
       
         
         output$confusionMatrix2<-DT::renderDataTable  ({
           #if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
               req(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)))
               mat=as.data.frame(t(confusion(clusterValues$comm2,clusterValues$comm1)$mat1))
               mat[,c(1,2)]<-mat[,c(2,1)]
               isolate({colnames(mat)<-c(input$clusterSelect2,input$clusterSelect1,"Freq")})
               mat
           
         })
         
         output$metrics2<-renderText({
           
              # if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
               req(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)))
               nmi=NMI(clusterValues$comm1,clusterValues$comm2)
               ari=adj.rand.index(clusterValues$comm1,clusterValues$comm2)
               paste0("NMI:",nmi," ARI:",ari)
           
         })
         
         
         output$interHeatmap<-renderPlot({
               #if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
               req(!(is.null(clusterValues$comm1) | is.null(clusterValues$comm2)))
               
               if (input$heatSelect==0) 
                  x=NULL
               else
                   x=as.numeric(input$heatSelect)
           
           
               heatmap.2(t(prop.table(table(clusterValues$comm1,clusterValues$comm2),x)), scale="none",trace = "none",
                         xlab="Clustering 1",ylab="Clustering 2",main="Probability Heatmap",density.info = "none",col=colorpanel(256,"blue","purple","red"),dendrogram = "none",Rowv = NA, Colv = NA)
               
                })
         
         
         output$downloadCluster <- downloadHandler(
               filename = function() {
                 paste(input$clusterSelect,".csv",sep="")
               },
               content = function(file) {
                 write.csv(clusterValues$member,file, row.names = FALSE)
               })    
         
         
         })

        
        
                
       
       
       
   
          
        
        
        
      
      
         
       
         
         
         
         
         
       
          
      


  
  
  

