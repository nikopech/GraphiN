library(shiny)
shinyUI(navbarPage("Navbar",
           tabPanel("Graph Creation",
        
        
                  fluidRow(
            
            
                    titlePanel("File Input"),
                    sidebarLayout(
                      sidebarPanel(
                        fileInput("file","Upload the file"), # fileinput() function is used to get the file upload contorl option
                        helpText("Default max. file size is 5MB"),
                        tags$hr(),
                        h5(helpText("Select the read.table parameters below")),
                        checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
                        br(),
                        radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ''),
                        br(),br(),br()
                        
                      ),
                       mainPanel(
                              
                                tabsetPanel(
                                    tabPanel("Graph",visNetworkOutput("network")),
                                    tabPanel("Data",numericInput("idInput","Insert the id",0),
                                                    dataTableOutput("dataset"))
                                )
                        # use below code if you want the tabset programming in the main panel. If so, then tabset will appear when the app loads for the first time.
                        #       tabsetPanel(tabPanel("Summary", verbatimTextOutput("sum")),
                        #                   tabPanel("Data", tableOutput("table")))
                      )
                      
                    )
                  ),
                 fluidRow(
                   br(),
                   selectInput("seqSelect","Select the column of the sequence",choices=list("IMGT.gapped.nt.sequences.V.D.J.REGION"="IMGT.gapped.nt.sequences.V.D.J.REGION",
                                                                                            "IMGT.gapped.AA.sequences.V.D.J.REGION"= "IMGT.gapped.AA.sequences.V.D.J.REGION",
                                                                                            "IMGT.gapped.nt.sequences.V.J.REGION"="IMGT.gapped.nt.sequences.V.J.REGION", 
                                                                                            "IMGT.gapped.AA.sequences.V.J.REGION"="IMGT.gapped.AA.sequences.V.J.REGION"),selected="IMGT.gapped.AA.sequences.V.D.J.REGION",width=350),
                   
                   selectInput("simSelect","Select the similarity metric",choices=list("OSA"="osa",
                                                                                       "Levenshtein"= "lv",
                                                                                       "Full DL"="dl", 
                                                                                       "Hamming"="hamming",
                                                                                       "LCS"="lcs",
                                                                                       "Q-Gram"="qgram",
                                                                                       "Cosine"="cosine",
                                                                                       "Jaccard"="jaccard",
                                                                                       "Jaro Winklar"="jw",
                                                                                       "Soundex"="soundex"),selected="osa",width=350),
                   
                   checkboxInput(inputId = "clusterId", label = "Unique sequence-clusterID combination", value = FALSE),
                   actionButton("simButton","Calculate distance matrix"),
                   tags$head(tags$script(src = "message-handler.js")),
                   sliderInput("slider", label = h3("Edge Threshold"), min = 0, max = 1, value = 0.6,step=0.001,width=500),
                   
                   br(),
                   actionButton("graphButton","Create Graph"),
                   br(),br(),
                   selectInput("componentSelect","Components",choices=list("no components"="")),
                   actionButton("componentButton","Select a component")
                 
                 
               )   
          
          
          ),
          tabPanel("GraphFiltering",
                   fluidRow(
                     column(2,
                        hr("Include"),
                        uiOutput("selectbox1"),
                        uiOutput("textbox1"),
                        actionButton("includeButton","Add to graph")
                        ),
                     column(4,
                        hr("Exclude"),
                        uiOutput("selectbox2"),
                        uiOutput("textbox2"),
                        actionButton("excludeButton","Delete from graph")
                        ),
                     column(6,
                       hr("Filters Added"),
                       tableOutput("filtertable")
                            )
                   ),
                   br(),br(),
                   fluidRow(
                      column(1,
                         actionButton("FilterButton","Filter")),
                      column(2,
                         actionButton("ResetButton","Reset")),
                      column(3,
                         checkboxInput("ReverseButton","Reverse Filters"))
                            ),
                   hr("Filtered Indexes"),
                   verbatimTextOutput("Indexes")
                    
                    ),
          
          
          tabPanel("MST",
                   
                   visNetworkOutput("mstnetwork"),
                   column(1,
                    plotOutput("mstLegend",width=500)),
                   column(2,
                    plotOutput("mstLegend2",width = 500)),  
                   selectInput("colormst","Select an attribute for background coloring",choices=c(Default="Default"),selected = "Default"),
                   selectInput("bordermst","Select an attribute for border coloring",choices=c(Default="Default"),selected = "Default"),
                   checkboxInput("clusterMST","Coloring according to clusters",value=FALSE),
                   actionButton("mstButton","Create MST"),
                   hr("Central Nodes"),
                   verbatimTextOutput("Cendroids"),
                   hr("Link nodes"),
                   verbatimTextOutput("Links")
            
            
          ),
          
          
          tabPanel("Centralities",
                   actionButton("centralButton","Calculate Centralities"),
                   br(),br(),
                   tabsetPanel(
                     
                          tabPanel("Centralities",DT::dataTableOutput("Centralities")),
                          tabPanel("Rankings",DT::dataTableOutput("Rankings")),
                          tabPanel("Summary",DT::dataTableOutput("Summary")),
                          tabPanel("Network",
                                   selectInput("centralSelect","Select a centrality type",list('Degree'="Degree.Centrality",
                                                                                                'Betweenness'="Shortest.Paths.Betweenness.Centrality",
                                                                                                'Closeness'="Average.Distance",
                                                                                                'Eigenvector'="eigenvector.centralities"),selected="Shortest.Paths.Betweenness.Centrality"),
                                                                                                
                                   visNetworkOutput("centralnetwork"),
                                   hr("Note:Triangle indicates the most 'central' code."))
                          )),
          
          
          tabPanel("Clustering",
                   tabsetPanel(
                            tabPanel("Graph with clusters",
                                  selectInput("clusterSelect","Select a clustering algorithm",list("Louvain"="louvain",
                                                                                     "Fast Greedy"="fast_greedy",
                                                                                     "Label Propagation"="label_propgation",
                                                                                     "Leading Eigenvalue"="leading_eigenvalue",
                                                                                     "Walktrap"="walktrap",
                                                                                     "Edge Betweeness"="edge_betweeness"),selected = "louvain"),
                                  fluidRow(
                                    column(4,
                                           visNetworkOutput("clusterNetwork")),
                                    column(8,
                                           plotOutput("intraHeatmap"))
                                  )),
                                  
                             tabPanel("Confusion Matrix",
                                  selectInput("clusterSelect1","Select a clustering algorithm",list("Louvain"="louvain",
                                                                                    "Fast Greedy"="fast_greedy",
                                                                                    "Label Propagation"="label_propgation",
                                                                                    "Leading Eigenvalue"="leading_eigenvalue",
                                                                                    "Walktrap"="walktrap",
                                                                                    "Edge Betweeness"="edge_betweeness"),selected = "louvain"),
                                  selectInput("clusterSelect2","Select a clustering algorithm",list("Louvain"="louvain",
                                                                                    "Fast Greedy"="fast_greedy",
                                                                                    "Label Propagation"="label_propgation",
                                                                                    "Leading Eigenvalue"="leading_eigenvalue",
                                                                                    "Walktrap"="walktrap",
                                                                                    "Edge Betweeness"="edge_betweeness"),selected = "louvain"),
                                
                                  fluidRow(
                                      column(4,
                                          tableOutput("confusionMatrix")),
                                      column(8,
                                          plotOutput("interHeatmap"))
                                          )
                                
                                     )
                )
           )
          
))
                
              
        
                      
                      
            
            
            
            
          
          
          
          
          
          
      
