library(shiny)
scatter_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    selectInput(
      "nb",
      "Number of points",
      choices = c(
        "genotype"   = 'genotype',
        "Organ"  = 'clusters',
        "expression"="expression"
      ),
      selected = 'clusters'
      
    ),
    conditionalPanel('output.show',selectizeInput("geneNameIn","gene:",choices=NULL)),
    plotlyOutput("p"),
    conditionalPanel(condition = "input.nb == 'clusters' |input.nb == 'genotype'",  plotlyOutput("bar"),
                     plotlyOutput("bar2")),
    conditionalPanel(condition = "input.nb == 'expression'",  plotlyOutput("dot")),
    
  )
}



scatter_server  <- function(id,color1,color4,genenames,merge14,testdata){
  moduleServer(
    id,
    function(input,output,session){
      
      ## gene selection
      updateSelectizeInput(session, 'geneNameIn', choices = genenames, server = TRUE,selected='egl-3')
      
      ## change color according to scatter or genotype
      output$p <- renderPlotly(
        plot_ly(source="scatter",testdata,x = ~UMAP_1, y = ~UMAP_2,
                mode = "markers",type="scattergl",opacity=0.5,marker=list(size=2,color=color4),
                customdata=~percent.mt,showlegend = T,
                hovertemplate=paste('<b>Mt</b>: %{customdata:.2f}%',
                                    '<br><b>Gene</b>: %{text}<extra></extra>'),
                text = ~paste( nFeature_RNA,'<br><b>RNA</b>:',nCount_RNA,'<br><b>Genotype</b>:',genotype,'<br><b>Cluster</b>:',seurat_clusters)
        ) %>% config(modeBarButtonsToRemove= c('toggleSpikelines','hoverCompareCartesian','toggleHover','hoverClosestCartesian','autoScale2d','select2d','lasso2d')) %>% 
          config(displaylogo=FALSE) %>% 
          config(scrollZoom= TRUE)  %>% 
          config(displayModeBar= TRUE) %>%
          event_register("plotly_hover") %>%
          layout(xaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 1'),
                 yaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 2'),
                 font=list(family='Arial',size=12),
                 dragmode='pan',title='Cells colored by Tissue')
      )
      
      #event_register(p, 'plotly_hover')
      observeEvent(input$nb, {
        if(input$nb=='clusters'){
          colori=c(color4)
          plotlyProxy("p", session) %>%
            plotlyProxyInvoke(
              "restyle", 
              list(marker=list(color=colori,size=2)) )}else if(input$nb=='genotype'){colori=c(color1)    
              plotlyProxy("p", session) %>%
                plotlyProxyInvoke(
                  "restyle", 
                  list(marker=list(color=colori,size=2)) )}
        
      })
      
      
      observeEvent(input$nb,{if(input$nb=='clusters'){titleout='Cells colored by Tissue'    
      plotlyProxy("p", session)  %>%
        plotlyProxyInvoke(
          "relayout", 
          list(title=titleout))}else if(input$nb=='genotype'){titleout='Cells colored by Genotype'
          plotlyProxy("p", session)  %>%
            plotlyProxyInvoke(
              "relayout", 
              list(title=titleout))}
        
      })
      
      ## interact with histogram
      output$bar <- renderPlotly(
        plot_ly(aa,x=~Var1,y=~Freq,type='bar',marker=list(color=myColors2,line=list(color='black',width=1.5))))
      
      output$bar2 <- renderPlotly(
        plot_ly(bb,x=~Var1,y=~Freq,type='bar',marker=list(color=myColors,line=list(color='black',width=1.5))))
      
      observeEvent(event_data("plotly_hover",source='scatter',session = session), {
        d <- event_data("plotly_hover",source='scatter',session=session)
        g1=testdata[d$pointNumber+1,]$genotype
        genotyepn=bb[bb$Var1==g1,]$X
        clusters=testdata[d$pointNumber+1,]$seurat_clusters
        clustern=aa[aa$Var1==clusters,]$X
        pattern_s_bar1<-rep(0.1,41)
        pattern_s_bar1[clustern]=1
        pattern_s_bar2<-rep(0.1,14)
        pattern_s_bar2[genotyepn]=1
        plotlyProxy("bar", session) %>%
          plotlyProxyInvoke(
            "restyle", 
            list(marker=list(opacity=pattern_s_bar1,color=myColors2,line=list(color='black',width=1.5)))
          )
        plotlyProxy("bar2", session) %>%
          plotlyProxyInvoke(
            "restyle", 
            list(marker=list(opacity=pattern_s_bar2,color=myColors,line=list(color='black',width=1.5)))
          )
      })  
      
      
      ## draw expression
      
      
      
      output$show <- reactive({input$nb=='expression'})
      
      outputOptions(output, 'show', suspendWhenHidden = FALSE)
      
      observeEvent(input$geneNameIn, {if(input$nb=='expression') {
        
        FetchData(object = merge14, vars = c(input$geneNameIn))->test1
        plotlyProxy("p", session) %>%
          plotlyProxyInvoke(
            "restyle", 
            list(marker=list(color=test1[,1],size=2)) )
        DotPlot(merge14,features = input$geneNameIn,split.by='genotype',cols = rep('red',14))$data->dotdata
        dotdata %>% separate(id, c("cluster", "genotype"),"_")->dotdata
        output$dot <- renderPlotly(plot_ly(dotdata,x=~cluster,y=~genotype,type='scatter',
                                           mode='markers',marker=list(size=~avg.exp,color=~pct.exp,line=list(color='black',width=2))))
        
      }})
      
      
    }
  )}








