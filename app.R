library(shiny)
library(dplyr)
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(tidyr)
library(shinymaterial)
library(Matrix)
library(ggpubr)
library(Seurat)


## data input

load("data/worm_seurat.RData")

source("modules/heatmap_main_8.R")


tissue_names_heatmap1 <- c("Intestine"
                           ,"Pharynx"
                           ,"Neuron"
                           ,"Hypodermis"
                           ,"Muscle"
                           ,"Spermatheca"
                           ,"Vulva & uterus"
                           ,"Germline")


tissue_names_heatmap2 <- c("Embroynic"
                  ,"Hypodermis"
                  ,"Germline"
                  ,"Spermatheca"
                  ,"Sperms"
                  ,"Muscle"
                  ,"Intestine"
                  ,"DTC/Excretory gland"
                  ,"Neuron"
                  ,"Gonadal sheath"
                  ,"Pharynx"
                  ,"Vulva & uterus"
                  ,"Uterine seam cell"
                  ,"Glia"
                  ,"Coelomocyte")

genotype_names_heatmap1 <- c("WT D1"
                    ,"WT D6"
                    ,"lipl-4 Tg D1"
                    ,"lipl-4 Tg D6"
                    ,"daf-2(If) D1"
                    ,"daf-2(If) D6"
                    ,"rsks-1(If) D1"
                    ,"rsks-1(If) D6")



cols<-brewer.pal(9, "RdYlBu")  
apa_umap_colors <- colorRampPalette(cols)




umap_colors_tissue <- c("#50AAAE" # Coelomocyte
  ,"#FDB668" # Embroynic
  ,"#5E4FA2" # DTC/Excretory gland
  ,"#FDD17E" # Germline
  ,"#358BBB" # Glia
  ,"#D4ED9B" # Gonadal sheath
  ,"#B5E1A1" # Hypodermis
  ,"#F0ADA0" # Intestine
  ,"#207F4C"# Muscle
  ,"#BB2148" # Neuron
  ,"#F8E48E" # Pharynx
  ,"#466DB0" # Spermatheca
  ,"#6EC5AC" # Sperms
  ,"#FFC0CB" # Uterine seam cell
  ,"#8abcd1" # Vulva & uterus
) 







hm_colors <- colorRampPalette(c("#FC4911", "#D83E0C", "#A92F08", 
                                "#741C03", "#3E0C00", "black", 
                                "#01233A", "#043250", "#0B4972", 
                                "#0F5B8D", "#1873B1"))
                                


geno_freq_colors <- c("#E5E5E5" # wt_D1
,"#BEBEBE" # wt_D6
,"#938F92" # wt_D12
,"#6B696B" # wt_D14
,"#6F84C1" # lipl4_D1
,"#0066CC" # lipl4_D6
,"#BCA15E" # daf2_D1
,"#CC9900" # daf2_D6
,"#9E5188" # rsks1_D1
,"#990066") # rsks1_D6


tissue_freq_colors <- c("#FDB668" # Embroynic
	,"#B5E1A1" # Hypodermis
	,"#FDD17E" # Germline
	,"#466DB0" # Spermatheca
	,"#6EC5AC" # Sperms
	,"#207F4C"# Muscle
	,"#F0ADA0" # Intestine
	,"#5E4FA2" # DTC/Excretory gland
	,"#BB2148" # Neuron
	,"#D4ED9B" # Gonadal sheath
	,"#F8E48E" # Pharynx
	,"#8abcd1" # Vulva & uterus
	,"#FFC0CB" # Uterine seam cell
	,"#358BBB" # Glia
	,"#50AAAE") # Coelomocyte






############## Functions pheatmap uses to carry out normalization ##############
scale_rowz <- function (x) 
{
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}


scale_matrix <- function (mat, scale) 
{
  if (scale == "row") {
    mat <- scale_rowz(mat)
    return(mat)
  } else if (scale == "column") {
    mat <- t(scale_rowz(t(mat)))
    return(mat)
  }
}
################################################################################

## data pre-processing

## UI
ui <-htmlTemplate(
  filename = "www/feb_20_template.html",
  expression = 
    
               
    fluidPage(
      fluidRow(column(5, selectInput('geneNameIn', 'Gene:',choices = 'ret-1'))),
      fluidRow(column(2),
               column(8,plotlyOutput("UMAP_expression")),
               column(2)),
      fluidRow(column(2),
               column(8,plotlyOutput("dot_expression")),
               column(2)),
      fluidRow(column(5, selectInput('hist', 'Number of Cells Grouped By:', choices = c('Tissue', 'Genotype')))),
      fluidRow(column(2),
               column(8,plotlyOutput("histograms")),
               column(2)

      
      )),
  
  APA = material_card(
    fluidPage(
      fluidRow(column(5, selectInput('geneNameIn_for_apa', 'Gene:', choices = 'ret-1'))),
      fluidRow(column(4, plotlyOutput("apa_N2D1_umap")),
               column(4, plotlyOutput("apa_N2D6_umap")),
               column(4, plotlyOutput("apa_N2D12.D14_umap"))),
      fluidRow(column(10), imageOutput("colorbar_output", height="50px")),
      fluidRow(column(5, selectInput("scale_by", "Scale heatmap by:", c("Genotype", "Tissue")))),
      fluidRow(column(8,plotlyOutput("heatmap1")),
               column(3,plotlyOutput("heatmap2")))
      ))
  )


################################################################################


server <- function(input, output,session) {
 
  updateSelectizeInput(session, 'geneNameIn', choices = genenames_seurat,server=TRUE,selected='ret-1')
  updateSelectizeInput(session, 'geneNameIn_for_apa', choices = genenames_APA, server = TRUE,selected='ret-1')
  

  
  observeEvent(input$geneNameIn, {

    FetchData(object=slimS,vars=c(input$geneNameIn))->expression_df
    
    output$UMAP_expression <- renderPlotly(
      
      
      plot_ly(source ='scatter',
              slimS[[]],x = ~UMAP_1, y = ~UMAP_2,
              mode = "markers",type="scattergl",opacity=0.5,marker=list(color=expression_df[,1],size=4),showlegend = T,
              
              hovertemplate=paste('<br><b>Tissue</b>: %{text}<extra></extra>'),
              text = ~paste( tissue,'<br><b>Genotype</b>:',genotype,'<br><b>Expression</b>:',round(expression_df[,1],digits = 2))
      ) %>% config(modeBarButtonsToRemove= c('toggleSpikelines','hoverCompareCartesian','toggleHover','hoverClosestCartesian','autoScale2d','select2d','lasso2d')) %>%
        config(displaylogo=FALSE) %>%
        config(scrollZoom= FALSE)  %>%

        layout(xaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 1'),
               yaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 2'),
               font=list(family='Arial',size=12),
               dragmode='pan',title='Cells Colored by Expression') %>%
        layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
               paper_bgcolor = "rgba(0, 0, 0, 0)")
    )
    
    
   
      DotPlot(slimS,features = input$geneNameIn,split.by = 'tissue',cols=rep('red',20))$data -> dotdata
      
      dotdata %>% separate(id,c("genotype","tissue"),"_") -> dotdata
      
      
      for (i in 1:nrow(dotdata)) {
        
        if (dotdata$genotype[i] == "LIPL4D1") {
          dotdata$genotype[i] = "lipl-4 Tg D1"
        } else if (dotdata$genotype[i] == "LIPL4D6") {
          dotdata$genotype[i] = "lipl-4 Tg D6"
        } else if (dotdata$genotype[i] == "N2D1") {
          dotdata$genotype[i] = "WT D1"
        } else if (dotdata$genotype[i] == "N2D6") {
          dotdata$genotype[i] = "WT D6"
        } else if (dotdata$genotype[i] == "N2D12") {
          dotdata$genotype[i] = "WT D12"
        } else if (dotdata$genotype[i] == "N2D14") {
          dotdata$genotype[i] = "WT D14"
        } else if (dotdata$genotype[i] == "DAF2D1") {
          dotdata$genotype[i] = "daf-2(If) D1"
        } else if (dotdata$genotype[i] == "DAF2D6") {
          dotdata$genotype[i] = "daf-2(If) D6"
        } else if (dotdata$genotype[i] == "RSKS1D1") {
          dotdata$genotype[i] = "rsks-1(If) D1"
        } else if (dotdata$genotype[i] == "RSKS1D6") {
          dotdata$genotype[i] = "rsks-1(If) D6"
        }
        
        if (dotdata$tissue[i] == "Embroynic.cells") {
          dotdata$tissue[i] = "Embroynic"
        } else if (dotdata$tissue[i] == "Excretory.gland.cell.and.dtc") {
          dotdata$tissue[i] = "DTC/Excretory gland"
        } else if (dotdata$tissue[i] == "Gonadal.sheath.cells") {
          dotdata$tissue[i] = "Gonadal sheath"
        } else if (dotdata$tissue[i] == "Vulva.uterus") {
          dotdata$tissue[i] = "Vulva & uterus"
        } else if (dotdata$tissue[i] == "Uterine.seam.cell") {
          dotdata$tissue[i] = "Uterine seam cell"
        } 
        
      }
      
      
      
      as.factor(dotdata$genotype)->dotdata$genotype
      

      
      output$dot_expression <-renderPlotly(plot_ly(dotdata,x=~tissue,y=~genotype,type='scatter',
                                                   mode='markers',marker=list(size=~avg.exp.scaled,color=~colors,line=list(color='black',width=2)),
                                                   hovertemplate=paste('<br><b>Average Expression</b>: %{text}<extra></extra>'),
                                                   text = ~paste( round(avg.exp, digits=2),'<br><b>Percent cells with expression</b>:',round(pct.exp, digits=2))
                                                   ) %>% layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",paper_bgcolor = "rgba(0, 0, 0, 0)")    %>%     config(displaylogo=FALSE) %>%
                                             config(scrollZoom= FALSE) 
                                           
                                           )
                                                   

  })
  
  
  ## UMAP_expression
  ##   histogram
  
  observeEvent(input$hist, {
    
    if (input$hist == 'Tissue') {
      
      output$histograms <- renderPlotly(
        plot_ly(tissue_freq,x=~TISSUE,y=~FREQUENCY,type='bar',
                marker=list(color=tissue_freq_colors,line=list(color='black',width=1.5))) %>%
          layout(xaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',
                              showline=TRUE,zeroline=F,title='Tissue'),
                 yaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',
                              showline=TRUE,zeroline=F,title='Cell Number',
                              type = "log", tickformat = ".1r"),
                 plot_bgcolor  = "rgba(0, 0, 0, 0)",
                 paper_bgcolor = "rgba(0, 0, 0, 0)",
                 dragmode=FALSE) %>%
          config(scrollZoom= FALSE) %>%
          config(displayModeBar= FALSE))
      
    } else if (input$hist == 'Genotype') {

output$histograms <- renderPlotly(
  plot_ly(genotype_freq,x=~GENOTYPE,y=~FREQUENCY,type='bar',
          marker=list(color=geno_freq_colors,line=list(color='black',width=1.5))) %>% 
    layout(xaxis = list(showgrid=FALSE, mirror=TRUE, ticks='outside',
                        showline=TRUE, zeroline=F, title='Genotype',
                        categoryorder = "array",
                        categoryarray = genotype_freq$GENOTYPE,
                        tickangle = -45), # specify the tick angle here
           yaxis = list(showgrid=FALSE, mirror=TRUE, ticks='outside',
                        showline=TRUE, zeroline=F, title='Cell Number'),
           plot_bgcolor = "rgba(0, 0, 0, 0)",
           paper_bgcolor = "rgba(0, 0, 0, 0)", 
           dragmode = FALSE) %>% 
    config(scrollZoom = FALSE) %>% 
    config(displayModeBar = FALSE))

    }


})




observeEvent(input$geneNameIn_for_apa, {
  
  gene = input$geneNameIn_for_apa
  
  n2d1_df <- N2D1_APA[,gene]
  n2d6_df <- N2D6_APA[,gene]
  n2d12.d14_df <- N2D12_D14_APA[,gene]
  
  
  cols<-brewer.pal(9, "RdYlBu")
  limz <- c(min(n2d1_df, n2d6_df, n2d12.d14_df), 
            max(n2d1_df, n2d6_df, n2d12.d14_df))
  
  n2d1_p <- ggplot(ident_data[which(ident_data$genotype == "N2D1"), ], aes(x=UMAP_1, y=UMAP_2, colour=n2d1_df , text = paste("APA preference: ", round(n2d1_df,2))  )) + 
    geom_point(size=rel(0.5), show.legend = FALSE) + 
    scale_colour_gradientn(colours = cols, na.value="grey95", limits = limz, name="", n.breaks = 10) +
    ggtitle("WT D1") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  n2d6_p <- ggplot(ident_data[which(ident_data$genotype == "N2D6"), ], aes(x=UMAP_1, y=ident_data[which(ident_data$genotype == "N2D6"),"UMAP_2"], colour=n2d6_df, text = paste("APA preference: ", round(n2d6_df,2)) )) + 
    geom_point(size=rel(0.5), show.legend = FALSE) + 
    scale_colour_gradientn(colours = cols, na.value="grey95", limits = limz) +
    ggtitle("WT D6") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  n2d12.d14_p <- ggplot(ident_data[which(ident_data$genotype == "N2D12" | ident_data$genotype == "N2D14"), ], aes(x=UMAP_1, y=UMAP_2, colour=n2d12.d14_df, text = paste("APA preference: ", round(n2d12.d14_df,2)) )) + 
    geom_point(size=rel(0.5), show.legend = FALSE) + 
    scale_colour_gradientn(colours = cols, na.value="grey95", limits = limz) +
    ggtitle("WT D12/D14") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  leg <- get_legend(n2d1_p + 
                    geom_point(show.legend = TRUE) +
                    theme(legend.direction ="horizontal",
                          legend.key.height = unit(0.5, 'cm'),
                          legend.key.width = unit(6, 'cm'),
                          legend.text = element_text(size=15))) 
                    
                                            
  
  
  output$colorbar_output <-  renderPlot(
    as_ggplot(leg) 
  )
  

  output$apa_N2D1_umap <- renderPlotly(
    
    ggplotly(n2d1_p,
            showlegend = T,tooltip = "text"
      ) %>% 
      config(modeBarButtonsToRemove = c(
        'toggleSpikelines',
        'hoverCompareCartesian',
        'toggleHover',
        'hoverClosestCartesian',
        'autoScale2d',
        'select2d',
        'lasso2d'
      )) %>% 
      config(displaylogo = FALSE) %>% 
      config(scrollZoom = FALSE) %>% 
      layout(
        xaxis = list(showgrid = FALSE, mirror = TRUE, ticks = 'outside', showline = TRUE, zeroline = F, title = 'UMAP 1'),
        yaxis = list(showgrid = FALSE, mirror = TRUE, ticks = 'outside', showline = TRUE, zeroline = F, title = 'UMAP 2'),
        font = list(family = 'Arial', size = 12),
        dragmode = 'pan') %>% 
      layout(
        plot_bgcolor = "rgba(0, 0, 0, 0)",
        paper_bgcolor = "rgba(0, 0, 0, 0)"
      )
  
)
  
  
  
  
  output$apa_N2D6_umap <- renderPlotly(
    ggplotly(n2d6_p,
             showlegend = T,tooltip = "text"

    ) %>% config(modeBarButtonsToRemove= c('toggleSpikelines','hoverCompareCartesian','toggleHover','hoverClosestCartesian','autoScale2d','select2d','lasso2d')) %>% 
      config(displaylogo=FALSE) %>% 
      config(scrollZoom= FALSE)  %>% 
      
      layout(xaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 1'),
             yaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 2'),
             font=list(family='Arial',size=12),
             dragmode='pan',title='WT D6') %>% layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
                                                                            paper_bgcolor = "rgba(0, 0, 0, 0)")
  )
  
  
  
  
  
 
  
  
  
  output$apa_N2D12.D14_umap <- renderPlotly(
    ggplotly(n2d12.d14_p,
             showlegend = T,tooltip = "text"
    ) %>% config(modeBarButtonsToRemove= c('toggleSpikelines','hoverCompareCartesian','toggleHover','hoverClosestCartesian','autoScale2d','select2d','lasso2d')) %>% 
      config(displaylogo=FALSE) %>% 
      config(scrollZoom= FALSE)  %>% 
       
      layout(xaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 1'),
             yaxis = list(showgrid=FALSE, mirror=TRUE,ticks='outside',showline=TRUE,zeroline=F,title='UMAP 2'),
             font=list(family='Arial',size=12),
             dragmode='pan',title='WT D12/D14') %>% 
      layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
             paper_bgcolor = "rgba(0, 0, 0, 0)")
  )
  
  
})




##### heatmaps #####
observeEvent(c(input$geneNameIn_for_apa,input$scale_by),{
  
  if (input$scale_by == "Genotype") {
    scale_input <- "column"
  } else if (input$scale_by == "Tissue") {
    scale_input <- "row"
  }
    
    df <- heatmap.main.8(input$geneNameIn_for_apa, avg_APA)
    
    output$heatmap1 <- renderPlotly({
      
      # Normalize the data by column
      APA <- scale_matrix(df, scale = scale_input)
      as.matrix(APA) -> APA
      
      
      plot_ly(z = ~APA, type = "heatmap", colors = hm_colors(150), 
              colorbar = list(title = "APA Score")) %>%
        layout(xaxis = list(title = list(text = "Genotypes", size = 20), 
                            ticktext = genotype_names_heatmap1, tickvals = 0:ncol(APA),
                            tickfont = list(size = 14)),
               yaxis = list(title = list(text = "Tissues", size = 20), 
                            ticktext = tissue_names_heatmap1, tickvals = 0:nrow(APA),
                            tickfont = list(size = 14)),
               font=list(family='Arial')) %>% 
        layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
               paper_bgcolor = "rgba(0, 0, 0, 0)") %>% 
        config(modeBarButtonsToRemove= c('toggleSpikelines',
                                         'hoverCompareCartesian',
                                         'toggleHover',
                                         'hoverClosestCartesian',
                                         'autoScale2d',
                                         'select2d',
                                         'lasso2d')) %>% 
        config(displaylogo=FALSE) %>% 
        config(scrollZoom= FALSE)

      
    })
    
    
  })



observeEvent(input$geneNameIn_for_apa, {
  
  
  df_heatmap2 <- avg_APA_N2D1[input$geneNameIn_for_apa,]
  rownames(df_heatmap2) <- "WT_D1"
  

  df_heatmap2 <- scale_matrix(df_heatmap2, scale = "row")
  as.matrix(t(df_heatmap2)) -> df_heatmap2
  
  
  output$heatmap2 <- renderPlotly({
    
    plot_ly(z = ~df_heatmap2, type = "heatmap", colors = hm_colors(150),
            width = 350, colorbar = list(title = "APA Score")) %>%
      layout(xaxis = list(title = list(text = "Genotype", size = 20),
                          ticktext = "WT D1", tickvals = 0:ncol(df_heatmap2),
                          tickfont = list(size = 14)),
             yaxis = list(title = list(text = "Tissues", size = 20),
                          ticktext = tissue_names_heatmap2, tickvals = 0:nrow(df_heatmap2),
                          tickfont = list(size = 14)),
             font=list(family='Arial')) %>% 
      layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
             paper_bgcolor = "rgba(0, 0, 0, 0)") %>% 
      config(modeBarButtonsToRemove= c('toggleSpikelines',
                                       'hoverCompareCartesian',
                                       'toggleHover',
                                       'hoverClosestCartesian',
                                       'autoScale2d',
                                       'select2d',
                                       'lasso2d')) %>% 
      config(displaylogo=FALSE) %>% 
      config(scrollZoom= FALSE) 
    
  })
  
  
  
})


}

################################################################################

shinyApp(ui = ui, server = server)








