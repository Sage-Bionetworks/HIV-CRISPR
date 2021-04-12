
library(shiny)
library(tidyverse)
library(DT)
library(plotly)
library(crosstalk)
library(PNWColors)

# Load cached data from Synapse
load("test_comp_data.RData")

ui <- fluidPage(
  
  titlePanel("HIV-CRISPR data viz demo: 20201230_Dragonite_Lat_CUL3_5A8v1WT"),
  
  tabsetPanel(
    tabPanel("Metadata",
             sidebarLayout(
               sidebarPanel(
                 "input TBD"
               ),
               mainPanel(
                 br(),
                 dataTableOutput("metadata")
               )
             )
             
    ),    
    tabPanel("Output data",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("showdata", "Select output files to show:",
                              choices = list("median_norm.gene_summary" = "median_norm",
                                             "control_norm.gene_summary" = "control_norm"),
                              selected = "median_norm"),
                 helpText("Use Ctrl/Cmd+click to select multiple rows. Search box accepts regular expressions (ie. use '|' for OR)")
               ),
               mainPanel(
                 br(),
                 dataTableOutput("comp_data_table")
               )
             )
             
    ),
    tabPanel("QC",
             sidebarLayout(
               sidebarPanel(
                 br(),
                 helpText("Uncheck the box to select specific sgRNAs 
                          (plot will only show with 3 or more selections)"),
                 checkboxInput("selectall", "Plot all sgRNAs", value = TRUE),
                 conditionalPanel(
                   condition = "input.selectall == false",
                   selectizeInput("pickgenes1", "Choose sgRNA(s) to plot:",
                                 choices = treatment_joined$sgRNA,
                                 multiple = TRUE,
                                 selected = c("CASP2_0", "CASP2_1", "CASP2_3"))
                 )
               ),
               mainPanel(
                 br(),
                 plotlyOutput("scatter_r2")
               )
             )
             
    ),
    tabPanel("Ranked gene summary",
             sidebarLayout(
               sidebarPanel(
                 sliderInput("topn", "Number of genes to plot:",
                             min = 1, max = 50, value = 20),
                 helpText("Rank by:"),
                 radioButtons("fillby", "Fill bars by:",
                              choices = list("gene in NTC list" = "ntc",
                                             "p <= 0.01" = "pval"),
                              selected = "ntc"),
                 helpText("FDR; score")
               ),
               mainPanel(
                 plotlyOutput("gene_bar")
               )
             )
             
    ),
    tabPanel("Individual sgRNAs",
             sidebarLayout(
               sidebarPanel(
                 selectizeInput("pickgenes2", "Choose genes to plot:",
                                choices = median_norm$id,
                                multiple = TRUE,
                                selected = "CASP2")
               ),
               mainPanel(
                 plotlyOutput("scatter_sgrna"),
                 p("Grey line indicates median of group")
               )
             )
             
    )
    
  )
  
)

server <- function(input, output) {
  
  # comparison metadata (from treatment replicate 1)
  output$metadata <- renderDataTable({
    metadata %>% 
      filter(id == "syn23702679") %>%
      select(-starts_with("ROW_")) %>% 
      pivot_longer(cols = everything(),
                   names_to = "field", values_to = "value",
                   values_transform = list(value = as.character)) %>% 
      datatable(caption = "Metadata for this comparison, pulled from treatment replicate 1")
  })
  

  # QC scatter plot with R2 value
  output$scatter_r2 <- renderPlotly({
    
    treatment_joined <- if (input$selectall == TRUE) {
      treatment_joined
    } else {
      treatment_joined %>% 
        filter(sgRNA %in% input$pickgenes1)
    }
    
    control_joined <- if (input$selectall == TRUE) {
      control_joined
    } else {
      control_joined %>% 
        filter(sgRNA %in% input$pickgenes1)
    }
    
    treatment_r2 <- round(cor.test(treatment_joined$Dragonite_20201201_1.fastq,
                                   treatment_joined$Dragonite_20201201_2.fastq,
                                   method = "pearson")$estimate ^ 2, 2)
    
      p1 <- treatment_joined %>%
      ggplot(aes(x = Dragonite_20201201_1.fastq, 
                 y = Dragonite_20201201_2.fastq,
                 dummy = sgRNA, group = 1)) +
      geom_point() +
      geom_smooth(method = "lm", 
                  formula = "y ~ x",
                  se = FALSE) +
      annotate("text", label = paste0("R^2 = ", treatment_r2), 
               x = min(treatment_joined$Dragonite_20201201_1.fastq), y
               = max(treatment_joined$Dragonite_20201201_2.fastq)) +
      labs(title = "Counts in treatment (left) and control (right) replicates")
    
    fig1 <- ggplotly(p1) %>% 
      style(textposition = "right")
    
    control_r2 <- round(cor.test(control_joined$Dragonite_20201201_17.fastq,
                                 control_joined$Dragonite_20201201_18.fastq,
                                 method = "pearson")$estimate ^ 2, 2)
    
    p2 <- control_joined %>% 
      ggplot(aes(x = Dragonite_20201201_17.fastq, 
                 y = Dragonite_20201201_18.fastq,
                 dummy = sgRNA, group = 1)) +
      geom_point() +
      geom_smooth(method = "lm", 
                  formula = "y ~ x",
                  se = FALSE) +
      annotate("text", label = paste0("R^2 = ", control_r2), 
               x = min(control_joined$Dragonite_20201201_17.fastq), 
               y = max(control_joined$Dragonite_20201201_18.fastq))
    
    fig2 <- ggplotly(p2) %>% 
      style(textposition = "right")
    
    subplot(fig1, fig2,
            titleX = TRUE,
            titleY = TRUE)
  })
  
  # choose dataset
  
  df <- reactive({
    if (input$showdata == "median_norm"){
      df <- median_norm
    } else {
      df <- control_norm
    }
  })
  
  
  # show data for single comparison
  output$comp_data_table <- renderDataTable({
    
    datatable(df(),
              options = list(search = list(regex = TRUE)),
              caption = paste0(input$showdata, ".gene_summary"))
    
  }) 
  
  # top 20 ranked genes
  output$gene_bar <- renderPlotly({
    
    df_top20 <- df() %>% 
      select(id, `neg|score`, `neg|rank`, 
             `pos|rank`, `pos|score`,
             `neg|p-value`, `pos|p-value`) %>%
      pivot_longer(cols = ends_with("rank"), names_to = "rank_type", values_to = "rank") %>% 
      pivot_longer(cols = ends_with("score"), names_to = "score_type", values_to = "score") %>% 
      pivot_longer(cols = ends_with("value"), names_to = "p_type", values_to = "p_value") %>% 
      arrange(rank) %>% 
      filter(rank <= input$topn) %>%
      filter(case_when(rank_type == "neg|rank" ~ score_type == "neg|score",
                       rank_type == "pos|rank" ~ score_type == "pos|score")) %>%
      filter(case_when(rank_type == "neg|rank" ~ p_type == "neg|p-value",
                       rank_type == "pos|rank" ~ p_type == "pos|p-value")) 
    
    
    p3 <- df_top20 %>% 
      ggplot(aes(y = fct_reorder(id, -rank), x = -log10(score),
                 fill = case_when(input$fillby == "ntc" ~ id %in% CUL3_synNTC_list,
                                  input$fillby == "pval" ~ p_value <= 0.01),
                 rank = rank, score = score, p = p_value)) +
      geom_col() +
      facet_wrap(~rank_type, scales = "free_y") +
      scale_fill_brewer(palette = "Dark2",
                        name = case_when(input$fillby == "ntc" ~ "gene in NTC list",
                                         input$fillby == "pval" ~ "p <= 0.01")) +
      ylab(NULL) +
      xlab("-log10 MAGeCK gene score") +
      ggtitle(paste0("Top ", input$topn, " highest-ranked genes in ", input$showdata, ".gene_summary"))
    
    ggplotly(p3, tooltip = c("rank", "score", "p"))
    
  })
  
  output$scatter_sgrna <- renderPlotly({
    
    genes_from_table <- df() %>% 
      filter(row_number() %in% input$comp_data_table_rows_selected) %>% 
      pull(id)
      
    
    p4 <- median_norm_sgRNA %>% 
      #filter(Gene %in% input$pickgenes2) %>% 
      filter(Gene %in% genes_from_table) %>% 
      ggplot(aes(x = Gene, y = score, color = Gene)) + 
      geom_jitter(aes(group = Gene, sgRNA = sgrna),
                  width = 0.1, alpha = 0.6) +
      stat_summary(aes(group = Gene),
                   fun = median, geom = "crossbar", 
                   width = 0.4, color = "darkgrey") +
      coord_flip() +
      scale_color_viridis_d() +
      theme(legend.position = "none") +
      ggtitle("Individual sgRNA scores for selected genes")
    
    ggplotly(p4, tooltip = c("sgRNA", "y"))
    
  })
  
}

shinyApp(ui = ui, server = server)
