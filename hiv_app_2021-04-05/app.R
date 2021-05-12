
library(shiny)
library(tidyverse)
library(tidytext)
library(DT)
library(plotly)
library(crosstalk)

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
                 helpText("Drag columns to reorder, 
                          or click 'Column visibility' to show/hide columns. 
                          Click to select multiple rows for plotting. 
                          Selected genes will be listed below. 
                          Search box accepts regular expressions (ie. use '|' for OR)"),
                 actionButton("clear", "Clear selections"),
                 br(),
                 br(),
                 wellPanel(
                   p("Selected genes: ", 
                     textOutput("selection_info1", inline = TRUE))
                 )
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
                 helpText("Uncheck the box to use gene selection from Data tab."),
                 checkboxInput("selectall_qc", "Plot all sgRNAs", value = TRUE),
                 br(),
                 conditionalPanel(
                   condition = "input.selectall_qc == 0",
                   wellPanel(
                     p("Selected genes: ", 
                       textOutput("selection_info2", inline = TRUE))
                   ))
               ),
               mainPanel(
                 br(),
                 plotlyOutput("scatter_r2")
               )
             )
             
    ),
    tabPanel("Individual sgRNAs",
             sidebarLayout(
               sidebarPanel(
                 helpText("Use the Data tab to select genes to plot. 
                          Vertical lines indicate median of group. 
                          (Error message means you haven't selected anything from the Data tab.)"),
                 br(),
                 wellPanel(
                   p("Selected genes: ", 
                     textOutput("selection_info3", inline = TRUE))
                 ),
                 radioButtons("sgrna_y", "Select variable to plot from sgRNA_summary output file:",
                              choices = list("Log2 fold change" = "LFC",
                                             "Score" = "score"),
                              selected = "LFC")
                 
               ),
               mainPanel(
                 plotlyOutput("scatter_sgrna")
               )
             )
             
    ),
    tabPanel("Ranked gene summary",
             sidebarLayout(
               sidebarPanel(
                 br(),
                 helpText("Uncheck the box to use gene selection from Data tab."),
                 checkboxInput("selectall_bar", "Rank all genes", value = TRUE),
                 br(),
                 # show selected genes if "rank all genes" is unchecked
                 conditionalPanel(
                   condition = "input.selectall_bar == 0",
                   wellPanel(
                     p("Plot selected genes: ", 
                       textOutput("selection_info4", inline = TRUE))
                   )
                 ),
                 # only show slider if "rank all genes" is checked
                 conditionalPanel(
                   condition = "input.selectall_bar == 1",
                   sliderInput("topn", "Number of genes to plot:",
                               min = 1, max = 50, value = 20)
                 ),
                 radioButtons("fillby", "Fill bars by:",
                              choices = list("gene in NTC list" = "ntc",
                                             "p <= 0.01" = "pval",
                                             "fdr <= 0.05" = "fdr"),
                              selected = "ntc")
               ),
               mainPanel(
                 conditionalPanel(
                   condition = "input.selectall_bar == 0",
                   plotlyOutput("genebar_selected")
                 ),
                 conditionalPanel(
                   condition = "input.selectall_bar == 1",
                   plotlyOutput("genebar_topn")
                 )
                 
               )
             )
             
    ),
    tabPanel("Dot Plot",
             sidebarLayout(
               sidebarPanel(
               ),
               mainPanel(
                 plotlyOutput("dotplot")
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
      datatable(caption = "Metadata for this comparison, pulled from treatment replicate 1",
                rownames = FALSE)
  })
  
  # choose dataset for Data tab
  df_gene <- reactive({
    if (input$showdata == "median_norm"){
      df_gene <- median_norm
    } else {
      df_gene <- control_norm
    }
  })
  
  # show data for single comparison
  output$comp_data_table <- renderDataTable({
    
    # bring in GeneCards links
    # don't replace gene id with link - there are a few differences between `MyList` and `Gene Symbol`
    df_gene_gc <- df_gene() %>% 
      left_join(select(CUL3_GO_GC, genecards, `Gene Symbol`), 
                by = c("id" = "Gene Symbol")) %>% 
      select(id, genecards, everything())
    
    datatable(df_gene_gc,
              rownames = FALSE,
              extensions = c("Buttons", "ColReorder"),
              options = list(search = list(regex = TRUE),
                             dom = "Bfrtip",
                             buttons = I("colvis"),
                             colReorder = list(realtime = FALSE)),
              caption = paste0(input$showdata, ".gene_summary"),
              escape = FALSE
    )
    
  }) 
  
  # set up proxy for "Clear Selections" button
  proxy <- dataTableProxy("comp_data_table")
  
  # action when "Clear Selections" button is clicked
  observeEvent(input$clear, {
    proxy %>% selectRows(NULL)
  })
  
  # pull vector of gene selections from datatable
  genes_from_table <- reactive({
    df_gene() %>%
      filter(row_number() %in% input$comp_data_table_rows_selected) %>%
      pull(id)
  })
  
  # print list of selected genes
  # do it multiple times bc can't reuse same output
  output$selection_info1 <- 
    output$selection_info2 <- 
    output$selection_info3 <-
    output$selection_info4 <-
    renderText({
      genes_from_table()
    },
    sep = ", ")
  
  # QC scatter plot with R2 value
  output$scatter_r2 <- renderPlotly({
    
    # plot all genes or a subset, based on input
    treatment_joined <- if (input$selectall_qc == TRUE) {
      treatment_joined
    } else {
      treatment_joined %>% 
        filter(str_detect(sgRNA, str_c(genes_from_table(), collapse = "|")))
    }
    
    control_joined <- if (input$selectall_qc == TRUE) {
      control_joined
    } else {
      control_joined %>% 
        filter(str_detect(sgRNA, str_c(genes_from_table(), collapse = "|")))
    }
    
    # calculate R2 and plot treatment replicates
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
                  se = FALSE) 
    
    fig1 <- ggplotly(p1) %>% 
      add_annotations(text = paste0("R^2 = ", treatment_r2),
                      x = 0, y = 1, xref = "x", yref = "paper", 
                      xanchor = "left", showarrow = FALSE)
    
    # calculate R2 and plot control replicates
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
                  se = FALSE) 
    
    fig2 <- ggplotly(p2) %>% 
      add_annotations(text = paste0("R^2 = ", control_r2), 
                      x = 0, y = 1, xref = "x", yref = "paper", 
                      xanchor = "left", showarrow = FALSE)
    
    # put the two plots together
    subplot(fig1, fig2,
            titleX = TRUE,
            titleY = TRUE) %>% 
      layout(title = "Counts in treatment (left) and control (right) replicates",
             margin = list(t = 50))
  })
  
  ## Two different plots for ranked gene summary: genebar_topn, genebar_selected
  
  # if ranking top N from all genes
  output$genebar_topn <- renderPlotly({
    
    df_top20 <- df_gene() %>%
      select(id, `neg|score`, `neg|rank`,
             `pos|rank`, `pos|score`,
             `neg|p-value`, `pos|p-value`,
             `neg|fdr`, `pos|fdr`) %>%
      pivot_longer(cols = ends_with("rank"), names_to = "rank_type", values_to = "rank") %>%
      pivot_longer(cols = ends_with("score"), names_to = "score_type", values_to = "score") %>%
      pivot_longer(cols = ends_with("value"), names_to = "p_type", values_to = "p_value") %>%
      pivot_longer(cols = ends_with("fdr"), names_to = "fdr_type", values_to = "fdr") %>%
      arrange(rank) %>%
      filter(rank <= input$topn) %>%
      filter(case_when(rank_type == "neg|rank" ~ score_type == "neg|score",
                       rank_type == "pos|rank" ~ score_type == "pos|score")) %>%
      filter(case_when(rank_type == "neg|rank" ~ p_type == "neg|p-value",
                       rank_type == "pos|rank" ~ p_type == "pos|p-value")) %>% 
      filter(case_when(rank_type == "neg|rank" ~ fdr_type == "neg|fdr",
                       rank_type == "pos|rank" ~ fdr_type == "pos|fdr"))
    
    p3 <- df_top20 %>%
      ggplot(aes(y = fct_reorder(id, -rank), x = -log10(score),
                 fill = case_when(input$fillby == "ntc" ~ id %in% CUL3_synNTC_list,
                                  input$fillby == "pval" ~ p_value <= 0.01,
                                  input$fillby == "fdr" ~ fdr <= 0.05),
                 rank = rank, score = score, p = p_value, fdr = fdr)) +
      geom_col() +
      facet_wrap(~rank_type, scales = "free_y") +
      scale_y_reordered() +
      scale_fill_brewer(palette = "Dark2",
                        name = case_when(input$fillby == "ntc" ~ "gene in NTC list",
                                         input$fillby == "pval" ~ "p <= 0.01",
                                         input$fillby == "fdr" ~ "FDR <= 0.05")) +
      ylab(NULL) +
      xlab("-log10 MAGeCK gene score") +
      ggtitle(paste0("Top ", input$topn, " highest-ranked genes in ", input$showdata, ".gene_summary"))
    
    ggplotly(p3, tooltip = c("rank", "score", "p", "fdr"))
    
  })  
  
  output$genebar_selected <- renderPlotly({
    
    df_gene_selected <-  df_gene() %>% 
      filter(id %in% genes_from_table()) %>% 
      select(id, `neg|score`, `neg|rank`, 
             `pos|rank`, `pos|score`,
             `neg|p-value`, `pos|p-value`,
             `neg|fdr`, `pos|fdr`) %>%
      pivot_longer(cols = ends_with("rank"), names_to = "rank_type", values_to = "rank") %>% 
      pivot_longer(cols = ends_with("score"), names_to = "score_type", values_to = "score") %>% 
      pivot_longer(cols = ends_with("value"), names_to = "p_type", values_to = "p_value") %>% 
      pivot_longer(cols = ends_with("fdr"), names_to = "fdr_type", values_to = "fdr") %>%
      filter(case_when(rank_type == "neg|rank" ~ score_type == "neg|score",
                       rank_type == "pos|rank" ~ score_type == "pos|score")) %>%
      filter(case_when(rank_type == "neg|rank" ~ p_type == "neg|p-value",
                       rank_type == "pos|rank" ~ p_type == "pos|p-value")) %>% 
      filter(case_when(rank_type == "neg|rank" ~ fdr_type == "neg|fdr",
                       rank_type == "pos|rank" ~ fdr_type == "pos|fdr"))
    
    
    # use {tidytext} to reorder within facets
    # create new id2 variable for reordering (adds "___"); use id for fill
    p4 <- df_gene_selected %>%
      mutate(rank_type = as.factor(rank_type),
             id2 = reorder_within(as.factor(id), -rank, rank_type)) %>% 
      ggplot(aes(y = id2, x = -log10(score),
                 fill = case_when(input$fillby == "ntc" ~ id %in% CUL3_synNTC_list,
                                  input$fillby == "pval" ~ p_value <= 0.01,
                                  input$fillby == "fdr" ~ fdr <= 0.05),
                 rank = rank, score = score, p = p_value, fdr = fdr)) +
      geom_col() +
      facet_wrap(~rank_type, scales = "free_y") +
      scale_y_reordered() +
      scale_fill_brewer(palette = "Dark2",
                        name = case_when(input$fillby == "ntc" ~ "gene in NTC list",
                                         input$fillby == "pval" ~ "p <= 0.01",
                                         input$fillby == "fdr" ~ "FDR <= 0.05")) +
      ylab(NULL) +
      xlab("-log10 MAGeCK gene score") +
      ggtitle(paste0("Selected genes from ", input$showdata, ".gene_summary, by pos or neg rank"))
    
    ggplotly(p4, tooltip = c("rank", "score", "p", "fdr"))
    
  })
  
  # choose dataset for sgRNA tab
  df_sgRNA <- reactive({
    if (input$showdata == "median_norm"){
      df_sgRNA <- median_norm_sgRNA
    } else {
      df_sgRNA <- control_norm_sgRNA
    }
  })
  
  # scatter plot of individual sgRNAs for each gene
  output$scatter_sgrna <- renderPlotly({
    
    p5 <- df_sgRNA() %>% 
      filter(Gene %in% genes_from_table()) %>% 
      ggplot(aes_string(x = "Gene", y = input$sgrna_y)) +
      # geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(group = Gene, sgRNA = sgrna, fill = Gene),
                  width = 0.1, alpha = 0.6, shape = 21) +
      stat_summary(aes(group = Gene),
                   fun = median, geom = "crossbar",
                   width = 0.5) +
      # coord_flip() +
      scale_fill_viridis_d() +
      theme(legend.position = "none") +
      xlab(NULL) +
      ggtitle("Individual sgRNA scores for selected genes")
    
    ggplotly(p5, tooltip = c("sgRNA", "y"))
    
  })
  
  # output$log2fc <- renderPlotly({
  #   
  #   p6 <- df_sgRNA() %>% 
  #     filter(Gene %in% genes_from_table()) %>% 
  #     ggplot(aes(x = Gene, y = LFC, group = Gene, sgrna = sgrna)) +
  #     geom_boxplot(outlier.shape = NA) +
  #     geom_jitter(aes(fill = Gene, group = Gene),
  #                 width = 0.1, alpha = 0.6, shape = 21) +
  #     scale_fill_viridis_d() +
  #     theme(legend.position = "none") +
  #     xlab(NULL) +
  #     ggtitle("Individual sgRNAs for selected genes")
  #   
  #   ggplotly(p6, tooltip = c("sgrna", "y"))
  #   
  # })
  
  output$dotplot <- renderPlotly({
    
    median_neg_ntc <- df_gene() %>% 
      filter(str_detect(id, "NTC")) %>% 
      pull(`neg|score`) %>% 
      median()
    
    median_pos_ntc <- df_gene() %>% 
      filter(str_detect(id, "NTC")) %>% 
      pull(`pos|score`) %>% 
      median()
    
    p7 <- df_gene() %>% 
      select(id, `neg|score`, `pos|score`) %>%
      mutate(neg_log_fold_change = log10(median_neg_ntc / `neg|score`),
             pos_log_fold_change = log10(median_pos_ntc / `pos|score`)) %>%
      # reverse the sign of pos column so it will be on the bottom of the plot
      mutate(pos_log_fold_change = -(pos_log_fold_change)) %>% 
      pivot_longer(cols = contains("log"),
                   names_to = "change_type",
                   values_to = "change_value") %>% 
      ggplot(aes(x = fct_shuffle(as.factor(id)),
                 y = change_value, 
                 fill = change_type, text = id)) +
      geom_point(shape = 21, alpha = 0.8, size = 2) +
      geom_point(data = . %>% filter(str_detect(id, "NTC")),
                 shape = 21, fill = "white", size = 2, 
                 show.legend = FALSE) +
      scale_fill_brewer(palette = "Dark2") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      xlab("Genes") +
      ylab("Log10 fold change from NTC median")
    
    ggplotly(p7, tooltip = "text")
    
  })
  
}

shinyApp(ui = ui, server = server)
