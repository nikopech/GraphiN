

# use the below options code if you wish to increase the file input limit, in this example file input limit is increased from 5MB to 9MB
# options(shiny.maxRequestSize = 9*1024^2)

shinyServer(function(input,output,session){
  
  values <- reactiveValues(flagEX=0,filterdf=data.frame(),ig=NULL,indexes=c(),forest=c())

  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
    fulldata <- eventReactive(input$file,{
        file1 <- input$file
        if(is.null(file1)){return()} 
        data=read.csv(file=file1$datapath, sep=input$sep, header = input$header)
        data=data[sample(max(dim(data)),100, replace = FALSE),]
        #sim=stringdistances(seq=as.character(data[,input$seqSelect]),algo=input$simSelect)
        values$indexes=1:nrow(data)
        data
       })
    
    
    sim<-eventReactive(input$simButton,{
      sim=stringdistances(seq=as.character(fulldata()[,input$seqSelect]),algo=input$simSelect)
      session$sendCustomMessage(type = 'testmessage',
                                message = 'Thank you for clicking')
      sim
      
      })
  
    
    observeEvent(input$graphButton, {
          
          templist=list(fulldata()[values$indexes,],sim()[values$indexes,values$indexes])
          
          
          
          
          
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
          
          V(ig)$value=V(ig)$freq_cluster_id^(1/3.5)*1
          #V(ig)$size=V(ig)$freq_cluster_id^(1/3.5)*5
          
          #V(ig)$color=visualiseGenes(data=V(ig)$dataName)$label
          
          V(ig)$color.background=colors[visualiseGenes(data=V(ig)$dataName)$label]
          V(ig)$color.border=V(ig)$color.background
          #shapes <- c("circle", "square", "csquare", "rectangle")
          #V(ig)$shape <- shapes[V(ig)$color]
          
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
          })
    
    
    observeEvent(input$componentButton,{
      if (is.null(values$ig)) {return()}
      values$ig=induced_subgraph(values$ig,which(values$forest==input$componentSelect),impl="auto")
      V(values$ig)$id=1:gorder(values$ig)
      values$forest=c()
      values$forest[1:gorder(values$ig)]=1
      updateSelectInput(session,"componentSelect",choices=list("1"=1),selected=1)
       })
  
        
        output$dataset<-DT::renderDataTable({
          temp=fulldata()
          if (is.null(temp)){return()}
          if (input$idInput!=0)
             {temp=temp[input$idInput==rownames(temp),]}
         
          temp
          
          })
    
    
    
        output$network <- renderVisNetwork({
          igtemp=values$ig
          if (is.null(igtemp)) {return()}
          data_vertices=as_data_frame(igtemp,what=c("vertices"))
          data_vertices=cbind(data_vertices,title=isolate(values$forest))
          data_edges=as_data_frame(igtemp,what=c("edges"))
          visNetwork(data_vertices, data_edges) %>%
            visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
            visOptions (nodesIdSelection = list("useLabels"=TRUE),selectedBy = list("variable"="title"),highlightNearest = TRUE)   %>%
            #visPhysics(solver="repulsion") %>%
            visInteraction(multiselect = TRUE)%>%
            visIgraphLayout(layout = "layout_with_fr",smooth=FALSE,physics = FALSE) %>%
            visLegend(addNodes = data.frame())
          #visClusteringByColor( colors) 
        })
          
          
          
       output$selectbox1<-renderUI({
         tempdata=fulldata()
         
         choice=list()
         for (l in colnames(tempdata))
            choice=c(choice,l)
         names(choice)=colnames(tempdata)
         selectInput("select1", label = h3("Select box"), choices = choice)
        })
       
       
       output$selectbox2<-renderUI({
         tempdata=fulldata()
         
         choice=list()
         for (l in colnames(tempdata))
           choice=c(choice,l)
         names(choice)=colnames(tempdata)
         selectInput("select2", label = h3("Select box"), choices = choice)
       })
       
         
      output$textbox1<-renderUI({
        tempdata=fulldata()
        selected=input$select1
        if (is.null(selected)) {return()}
        if (is.numeric(tempdata[,selected]))
           {sliderInput("slider1", label = h3("Slider Range"), min = min(tempdata[,selected]), max = max(tempdata[,selected]), value = c(min(tempdata[,selected]),max(tempdata[,selected]) ))}
        else
            {textInput("text1","Insert Text")}
        
        
        
        
        
      }) 
      output$textbox2<-renderUI({
        tempdata=fulldata()
        selected=input$select2
        if (is.null(selected)) {return()}
        if (is.numeric(tempdata[,selected]))
        {sliderInput("slider2", label = h3("Slider Range"), min = min(tempdata[,selected]), max = max(tempdata[,selected]), value = c(min(tempdata[,selected]),max(tempdata[,selected]) ))}
        else
        {textInput("text2","Insert Text")}
      })
        
 
      
      

      
      xxchange <- reactive({
        paste(input$excludeButton ,input$includeButton)
      })
     
      
      
      observeEvent(xxchange(),
                              {
                                if(input$excludeButton[[1]]==0 & input$includeButton[[1]]==0) {return()}
                                tempdata=fulldata()  
                                if (isolate(values$flagEX)!=input$excludeButton[[1]]){
                                    selected=input$select2
                                    
                                    if (is.numeric(tempdata[,selected])){
                                        x=data.frame("Columns"=selected,"Keys"=paste0(input$slider2[1],",",input$slider2[2]),"I/E"="E")}
                                    else{   
                                        x=data.frame("Columns"=selected,"Keys"=input$text2,"I/E"="E")}
                                    
                                    filterg<<-rbind(filterg,x)
                                    values$flagEX=input$excludeButton[[1]]
                                    
                                }
                                else 
                                  {
                                    selected=input$select1
                                    if (is.numeric(tempdata[,selected])){ 
                                        x=data.frame("Columns"=selected,"Keys"=paste0(input$slider1[1],",",input$slider1[2]),"I/E"="I")}
                                    else{
                                        x=data.frame("Columns"=selected,"Keys"=input$text1,"I/E"="I")}
                                    
                                    filterg<<-rbind(filterg,x)
                                    }
                                 values$filterdf<-filterg
                                
                              })
      
      
      observeEvent(input$ResetButton,{
                    filterg<<-data.frame(Columns=character(),Keys=character(),"I/E"=character())
                    values$filterdf<-filterg
                    })
      
      
      
      observeEvent(input$FilterButton,{
                  tempdata=isolate(fulldata())
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
      
      
      
      mstValues<-reactiveValues(ig=NULL,centroids=NULL,keyvertices=NULL,clusters=NULL,flag2=NULL,legenddf=NULL,legenddf2=NULL)
      
      
      observeEvent(input$mstButton,{
          
          if (is.null(values$ig)) {return()}
          mstig=mst(values$ig,algorithm='prim')
          a=mstClustering(values$ig)
          
         
          if ((input$clusterMST) & !is.null(clusterValues$member) )
          {V(mstig)$color.background=colors[clusterValues$member]
          mstValues$legenddf=data.frame("label"=unique(clusterValues$member),"color"=(colors[unique(clusterValues$member)]))}
            
          else
            {  
            if (input$colormst=="Default")
                {V(mstig)$color.background="black"  #colors[a$clusters]
                 mstValues$legenddf=NULL}#data.frame("label"=unique(a$clusters),"color"=unique(colors[a$clusters]))
            else
                {temp=visualiseGenes(get.vertex.attribute(mstig,input$colormst))
                 V(mstig)$color.background=colors[temp$label]
                 mstValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[unique(temp$label)])}}
          
          if (input$bordermst=="Default")
            {V(mstig)$color.border="black" #colors[a$clusters]
             mstValues$legenddf2=NULL} #data.frame("label"=unique(a$clusters),"color"=unique(colors[a$clusters]))}
          else
            {temp=visualiseGenes(get.vertex.attribute(mstig,input$bordermst))
             V(mstig)$color.border=colors[temp$label]
             mstValues$legenddf2=data.frame("label"=unique(temp$name),"color"=colors[unique(temp$label)])}
          
          V(mstig)$value=V(values$ig)$freq_cluster_id^(1/3.5)*1
          
          m=(1- E(mstig)$weight-min(1-E(mstig)$weight))/(max(1-E(mstig)$weight)-min(1-E(mstig)$weight))
          E(mstig)$width <- m*2.5
          n=c()
          n[gorder(mstig)+1]=1
          n[a]="triangle"
          V(mstig)$shape=n[1:gorder(mstig)]
          V(mstig)$id=1:gorder(mstig)
          mstValues$ig=mstig
          mstValues$centroids=V(mstig)$label[a] #a$centroids]
          #mstValues$keyvertices=V(mstig)$label[unique(unlist(a$keyvertices))]
          #mstValues$clusters=a$clusters
          mstValues$flag2=NULL
        
      },ignoreInit = TRUE)
      
      
      output$mstnetwork<-renderVisNetwork({
        
        if (is.null(mstValues$ig)) {return()}
        data_vertices=as_data_frame(mstValues$ig,what=c("vertices"))
        data_edges=as_data_frame(mstValues$ig,what=c("edges"))
        coords=layout_with_stress(mstValues$ig)
        
        mstValues$flag2=c(1)
        visNetwork(data_vertices, data_edges,groups=sample(c("A","B"),25,replace = TRUE)) %>%
          visNodes(borderWidth = 2) %>%
          visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
          visOptions (nodesIdSelection = TRUE,highlightNearest = TRUE)   %>%
          ##visPhysics(solver="forceAtlas2Based",forceAtlas2Based = list(gravitationalConstant=-2000,centralGravity=0.02,springLength=40,springConstant=0.4,damping=1,avoidOverlap=1)) %>%
          visInteraction(multiselect = TRUE)%>%
          visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) 
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
        updateSelectInput(session,"colormst", label = "Select a Column for background", choices = choice)
        updateSelectInput(session,"bordermst", label = "Select a column for border", choices = choice)
        updateSelectInput(session,"clusterSelect2",label="Select a clustering algorithm or a column",choices=c(choice[-1],list("Louvain"="louvain",
                                                                                                                           "Fast Greedy"="fast_greedy",
                                                                                                                           "Label Propagation"="label_propgation",
                                                                                                                           "Leading Eigenvalue"="leading_eigenvalue",
                                                                                                                           "Walktrap"="walktrap",
                                                                                                                           "Edge Betweeness"="edge_betweeness")),selected = "louvain")
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
                  mstValues$legenddf=data.frame("label"=unique(temp$name),"color"=colors[unique(temp$label)])}}
          
          if (y=="Default")
            {bordcolor="black" #colors[mstValues$clusters]
             mstValues$legenddf2=NULL}#data.frame("label"=unique(mstValues$clusters),"color"=unique(colors[mstValues$clusters]))}
          else  
            {temp=visualiseGenes(tempdata[,input$bordermst])
             bordcolor=colors[temp$label]
             mstValues$legenddf2=data.frame("label"=unique(temp$name),"color"=colors[unique(temp$label)])}
            
          visNetworkProxy("mstnetwork") %>%
              visUpdateNodes(data.frame(id=1:nrow(tempdata),color.background=bgcolor,color.border=bordcolor))
          
        })
        
        centralValues=reactiveValues(centralities=data.frame(),rankings=data.frame,summary=data.frame())
        
        observeEvent(input$centralButton,{
            if (is.null(values$ig)){return()}
            igtemp=giant_component_extract(values$ig)[[1]]
            prop=proper_centralities(igtemp)
            central=calculate_centralities(igtemp,include=prop[c(-2,-27,-5,-6,-8,-21,-43,-38:-39,-31:-33)])
            central=central[lengths(central)!=0]
            x=as.data.frame(central,col.names = names(central),row.names=V(igtemp)$label)
           
            central_ranked=as.data.frame(apply(x,2,sort))
            colnames(central_ranked)=colnames(x)
          
            b=c("Average.Distance","Local.Bridging.Centrality","Wiener.Index.Centrality")
            
            #revert the order of indices 
            central_ranked[,-match(b,colnames(central_ranked))]=central_ranked[nrow(central_ranked):1,-match(b,colnames(central_ranked))]
          
            positions=c()
            for (i in 1:length(x))
            {positions_first=match(x[,i],central_ranked[,i])
            positions_last=match(x[,i],central_ranked[nrow(x):1,i])
            positions=cbind(positions,(positions_first-positions_last+nrow(x)+1)/2)
            }
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
          
          
         })
      
         
         output$Summary<-DT::renderDataTable({
           
             if (is.null(centralValues$summary)) {return()}
             centralValues$summary
           
          })
         
         output$centralnetwork<-renderVisNetwork({
           
              igtemp=giant_component_extract(isolate(values$ig))[[1]]
              if (is.null(igtemp) | is.null(centralValues$centralities)) {return()}
              
              
              bc <- centralValues$centralities[,input$centralSelect]
              if (input$centralSelect=="Average.Distance") {bc=1/bc}
              p1 <- ggraph(igtemp,layout = "centrality", cent = bc,iter = 500, tol = 1e-04) +
                draw_circle(use = "cent") +
                annotate_circle(bc,format="",pos="bottom") +
                geom_edge_link(edge_color="black",edge_width=0.3)+
                geom_node_point(aes(fill=as.factor(Faction)),size=2,shape=21)+
                scale_fill_manual(values=c("#8B2323", "#EEAD0E"))+
                theme_graph()+
                theme(legend.position = "none")+
                coord_fixed()+
                labs(title="betweenness centrality")
              
              coords=cbind(x=p1$data$x,y=p1$data$y)
              
              i=c()
              i[gorder(igtemp)+1]=""
              i[coords[,1]==0 & coords[,2]==0]="triangle"
              
              data_vertices=cbind(as_data_frame(igtemp,what=c("vertices")),shape=i[1:gorder(igtemp)])
              data_vertices$id=rownames(data_vertices)
              data_edges=as_data_frame(igtemp,what=c("edges"))
              
              visNetwork(data_vertices, data_edges) %>%
                visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
                visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
                #visPhysics(solver="repulsion") %>%
                visInteraction(multiselect = TRUE)%>%
                visIgraphLayout(layout = "layout.norm", layoutMatrix = coords,smooth=FALSE,physics = FALSE) 
          })
         
         
         
         clusterValues=reactiveValues(member=NULL,comm1=NULL,comm2=NULL)
         
         
         output$clusterNetwork<-renderVisNetwork({
           igtemp=values$ig
           if (is.null(igtemp)){return()}
           model=communities_general(igtemp,algorithm =input$clusterSelect)
           clusterValues$member=model$membership
           V(igtemp)$grp=model$membership
           bb <- layout_as_backbone(igtemp,keep=0.4)
           V(igtemp)$color=colors[model$membership]
           
           data_vertices=cbind(as_data_frame(igtemp,what=c("vertices")))
           data_vertices$id=rownames(data_vertices)
           data_edges=as_data_frame(igtemp,what=c("edges"))
           
           
           visNetwork(data_vertices, data_edges) %>%
             visEdges(color=list(color="black",highlight="red"),labelHighlightBold=FALSE)%>%
             visOptions (nodesIdSelection = list("useLabels"=TRUE),highlightNearest = TRUE)   %>%
             #visPhysics(solver="repulsion") %>%
             visInteraction(multiselect = TRUE)%>%
             visIgraphLayout(layout = "layout.norm", layoutMatrix = cbind(x=bb$xy[,1],y=bb$xy[,2]),smooth=FALSE,physics = FALSE) 
           })
         
         
         
         output$intraHeatmap<-renderPlot({
              dis=as.matrix(isolate(values$ig[]))
              maxt=max(apply(dis,1,max))
              dis[dis==0]=maxt+0.1
              diag(dis)=0
              temp=order(clusterValues$member)
              dis=dis[temp,temp]
              heatmap(dis, Rowv = NA, Colv = NA, col = heat.colors(256), revC = TRUE)
              #model$membership  
           
           
           
           
           
         })
         
         
         
         observe({
            igtemp=values$ig
            if (is.null(igtemp)){return()}
            clusterValues$comm1=communities_general(igtemp,algorithm =input$clusterSelect1)$membership
            if (input$clusterSelect2 %in% c("louvain","fast_greedy","label_propgation","leading_eigenvalue","walktrap","edge_betweeness"))
                {clusterValues$comm2=communities_general(igtemp,algorithm =input$clusterSelect2)$membership}
            else {
                 clusterValues$comm2=get.vertex.attribute(igtemp,input$clusterSelect2)
              
              
              }
            })
         
         output$confusionMatrix<-renderTable({
            if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
            confusion(clusterValues$comm1,clusterValues$comm2)$mat1
            })
         
         
         output$interHeatmap<-renderPlot({
           if (is.null(clusterValues$comm1) | is.null(clusterValues$comm2)) {return()}
           heatmap( t(prop.table(table(clusterValues$comm1,clusterValues$comm2))), Rowv = NA, Colv = NA, col = heat.colors(256)[256:1], revC = TRUE)
           
            })
         
         
         # output$confusionMatrix<-renderTable({
         #   igtemp=values$ig
         #   if (is.null(igtemp)){return()}
         #   model1=communities_general(igtemp,algorithm =input$clusterSelect1)
         #   model2=communities_general(igtemp,algorithm =input$clusterSelect2)
         #   as.matrix(confusion(model1$membership,model2$membership)$mat1)
         #   
         #   
         # })
         
      
          
      
      })

    
        
       
       
       
   
          
        
        
        
      
      
         
       
         
         
         
         
         
       
          
      


  
  
  

