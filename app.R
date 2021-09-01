
## HIV-CRISPR DATA VIZ SHINY APP

library(shiny)
library(synapser)
library(synapserutils)
library(tidyverse)
library(tidytext)
library(DT)
library(plotly)
library(waiter)


##========================================================================================

ui <- fluidPage(
  
  # Synapse stuff from template
  tags$head(
    singleton(
      includeScript("www/readCookie.js")
    )
  ),
  
  # TITLE PANEL
  titlePanel("HIV-CRISPR Screen Data Visualization"),
  
  ##========================================================================================
  
  # SINGLE SCREEN DATA TAB
  # Display metadata for all screens and output data for 1 selected screen
  # Select a screen and genes of interest
  tabsetPanel(
    tabPanel("Single screen data",
             sidebarLayout(
               sidebarPanel(
                 helpText("Note: May take a moment to load",
                          br(), br(),
                          "Top table shows metadata for all available screens. 
                        Select a screen to visualize.",
                        br(), br(),
                        "Bottom table shows gene summary output data
                        (median or control normalized) for selected screen. 
                        Select one or more genes to visualize.",
                        br(), br(),
                        "Drag columns to reorder, 
                        or click 'Column visibility' to show/hide columns. 
                        Click to select multiple rows for plotting. 
                        Selected genes will be listed below. 
                        Search box accepts regular expressions (ie. use '|' for OR).",
                        br(), br()),
                 radioButtons("showdata", "Select output type to show:",
                              choices = list("median_norm.gene_summary" = "median_norm",
                                             "control_norm.gene_summary" = "control_norm"),
                              selected = "median_norm"),
                 actionButton("clear", "Clear selections"),
                 br(), br(),
                 wellPanel(
                   p("Selected genes: ", 
                     textOutput("selection_info1", inline = TRUE))
                 )
               ),
               mainPanel(
                 br(),
                 h4("Metadata for all available screens"),
                 dataTableOutput("all_screens"),
                 br(), hr(), br(),
                 h4("Output data for selected screen"),
                 dataTableOutput("comp_data_table")
               )
             )
             
    ),
    
    ##========================================================================================    
    
    # QC TAB
    # Plot treatment & control count files for selected screen
    # Toggle all genes or only selected genes
    tabPanel("QC",
             sidebarLayout(
               sidebarPanel(
                 br(),
                 helpText("Note: May take a moment to load after selecting a new dataset",
                          br(), br(),
                          "Uncheck the box to use gene selection from Data tab."),
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
                 h4("Correlation between treatment (left) and control (right) 
                    replicate count files"),
                 plotlyOutput("scatter_r2")
               )
             )
             
    ),
    
    ##========================================================================================
    
    # INDIVIDUAL sgRNAs TAB
    # Plot sgRNA_summary data for selected genes
    tabPanel("Individual sgRNAs",
             sidebarLayout(
               sidebarPanel(
                 helpText("Use the Data tab to select genes to plot. 
                          Horizontal bar indicates median of group."),
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
                 br(),
                 h4("Individual sgRNA results for selected genes"),
                 plotlyOutput("scatter_sgrna")
               )
             )
             
    ),
    
    ##========================================================================================
    
    # RANKED GENE SUMMARY TAB
    # Plot neg & pos gene scores, ordered by rank
    # Toggle all genes or selected genes
    # Choose filters/fills
    tabPanel("Ranked gene summary",
             sidebarLayout(
               sidebarPanel(
                 br(),
                 helpText("Uncheck the box to use gene selection from Data tab. 
                          (Error message means you haven't selected anything 
                          from the Data tab.)"),
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
                 h4("Negative (left) and positive (right) gene scores"),
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
    
    ##========================================================================================
    
    # DOT PLOT TAB
    # Plot all genes from selected screen - NTCs in white
    tabPanel("Dot plot",
             sidebarLayout(
               sidebarPanel(
                 helpText("Note: This plot only works with median_norm.gene_summary.txt. 
                          Log fold change for positive scores 
                          is arbitrarily shown below 0 for plotting purposes. 
                          NTCs are shown in white.")
               ),
               mainPanel(
                 h4("Log fold change from NTC median for all genes in selected screen"),
                 plotlyOutput("dotplot")
               )
             )
             
    ),
    
    ##========================================================================================
    
    # COMPARE 2 SCREENS TAB
    # Select 2 screens
    # Plot comparison of neg & pos scores
    # Click a point to show gene data across all screens
    tabPanel("Compare 2 Screens",
             sidebarLayout(
               sidebarPanel(
                 helpText("Top table shows metadata for all available screens. 
                        Select 2 screens to compare. 
                        (Error means you've only selected 1 screen.)"),
                 helpText("NTCs will only be shown when 
                          comparing two screens with the same library."),
                 br(),
                 p(strong("Selected screens: "),
                   br(),
                   htmlOutput("selected_2_screens", inline = TRUE)),
                 br(),
                 radioButtons("norm_type", "Select output files to use:",
                              choices = list("median_norm.gene_summary" = "median_norm",
                                             "control_norm.gene_summary" = "control_norm"),
                              selected = "median_norm"),
               ),
               mainPanel(
                 br(),
                 h4("Metadata for all available screens"),
                 br(),
                 dataTableOutput("all_screens_2"),
                 hr(),
                 h4("Comparison of negative (left) and positive (right) 
                    -log10(score) for 2 selected screens"),
                 plotlyOutput("compare_plot"),
                 br(), hr(),
                 h4("Data for selected gene across all screens"),
                 br(),
                 dataTableOutput("gene_all_screens_table")
               )
             )
             
    )
  ),
  
  ##========================================================================================
  
  # more Synapse template stuff
  # make sure this is outside the tabSetPanel
  use_waiter(),
  waiter_show_on_load(
    html = tagList(
      img(src = "loading.gif"),
      h4("Retrieving Synapse information...")
    ),
    color = "#424874"
  )
  
)


##========================================================================================

server <- function(input, output, session) {
  
  # Synapse stuff from template
  # shows fancy loading screen
  # checks if you are logged in to Synapse
  session$sendCustomMessage(type="readCookie", message=list())
  
  observeEvent(input$cookie, {
    # If there's no session token, prompt user to log in
    if (input$cookie == "unauthorized") {
      waiter_update(
        html = tagList(
          img(src = "synapse_logo.png", height = "120px"),
          h3("Looks like you're not logged in!"),
          span("Please ", a("login", 
                            href = "https://www.synapse.org/#!LoginPlace:0", 
                            target = "_blank"),
               " to Synapse, then refresh this page.")
        )
      )
    } else {
      # login and update session; otherwise, notify user to log in to Synapse first
      tryCatch({
        synLogin(sessionToken = input$cookie, rememberMe = FALSE)
        
        ### update waiter loading screen once login successful
        waiter_update(
          html = tagList(
            img(src = "synapse_logo.png", height = "120px"),
            h3(sprintf("Welcome, %s!", synGetUserProfile()$userName))
          )
        )
        Sys.sleep(2)
        waiter_hide()
      }, error = function(err) {
        Sys.sleep(2)
        waiter_update(
          html = tagList(
            img(src = "synapse_logo.png", height = "120px"),
            h3("Login error"),
            span(
              "There was an error with the login process. Please refresh your Synapse 
              session by logging out of and back in to",
              a("Synapse", href = "https://www.synapse.org/", target = "_blank"),
              ", then refresh this page."
            )
          )
        )
        
      })
      
      # Any shiny app functionality that uses Synapse should be within the
      # input$cookie observer (so basically everything)
      output$title <- renderUI({
        titlePanel(sprintf("Welcome, %s", synGetUserProfile()$userName))
        
      })
      
      ##========================================================================================      
      
      # SINGLE SCREEN DATA TAB
      
      # Load data from Synapse
      source("get_synapse_data.R")
      
      # Show metadata for all screens with ACCEPTED yaml file
      output$all_screens <- renderDataTable({
        metadata %>% 
          datatable(rownames = FALSE,
                    filter = "top",
                    extensions = c("Buttons", "ColReorder"),
                    options = list(search = list(regex = TRUE),
                                   dom = "Bfrtip",
                                   buttons = I("colvis"),
                                   colReorder = list(realtime = FALSE)),
                    escape = FALSE,
                    selection = list(mode = "single")) 
      })
      
      # get synID for selected screen
      sample_sheet_id <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_rows_selected) %>%  
          pull(configId)
      })
      
      # get name of selected screen
      selected_screen <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_rows_selected) %>%  
          pull(name)
      })
      
      # pull library name for eventual use in identifying NTCs
      library_name <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_rows_selected) %>%  
          pull(LibraryName)        
      })
      
      # choose dataset for single screen plots
      df_gene <- reactive({
        if (input$showdata == "median_norm"){
          df_gene <- get(paste0(selected_screen(), "_median_norm.gene_summary.txt"))
        } else {
          df_gene <- get(paste0(selected_screen(), "_control_norm.gene_summary.txt"))
        }
      })
      
      # show data for single screen
      output$comp_data_table <- renderDataTable({
        
        # suppress error message that shows when nothing is selected 
        req(input$all_screens_rows_selected)
        
        # bring in GeneCards links
        # don't replace gene id with link - `MyList` may be different from `Gene Symbol`
        ### Update when I get Metascape data for the rest of the libraries
        df_gene_gc <- df_gene() %>% 
          left_join(select(CUL3_GO_GC, genecards, `Gene Symbol`), 
                    by = c("id" = "Gene Symbol")) %>% 
          select(id, genecards, everything())
        
        datatable(df_gene_gc,
                  rownames = FALSE,
                  filter = "top",
                  extensions = c("Buttons", "ColReorder"),
                  options = list(search = list(regex = TRUE),
                                 dom = "Bfrtip",
                                 buttons = I("colvis"),
                                 colReorder = list(realtime = FALSE)),
                  caption = paste0(input$showdata, ".gene_summary for ", 
                                   selected_screen()),
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
      
      # print list of selected genes in sidebar
      # do it multiple times bc you can't reuse same output
      output$selection_info1 <- 
        output$selection_info2 <- 
        output$selection_info3 <-
        output$selection_info4 <-
        renderText({
          
          # suppress error message if nothing selected
          req(input$comp_data_table_rows_selected)
          
          genes_from_table()
        },
        sep = ", ")
      
      ##========================================================================================      
      
      # QC TAB
      
      # QC scatter plot of count files with R2 values
      output$scatter_r2 <- renderPlotly({
        
        # suppress error message if nothing selected
        req(input$all_screens_rows_selected)
        
        # get count files for selected screen
        get_counts(sample_sheet_id())
        
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
        treatment_r2 <- round(cor.test(treatment_joined[[2]],
                                       treatment_joined[[3]],
                                       method = "pearson")$estimate ^ 2, 2)
        
        p1 <- treatment_joined %>%
          ggplot(aes(x = treatment_joined[[2]], 
                     y = treatment_joined[[3]],
                     dummy = sgRNA, group = 1)) +
          geom_point() +
          geom_smooth(method = "lm", 
                      formula = "y ~ x",
                      se = FALSE) +
          xlab(colnames(treatment_joined[2])) +
          ylab(colnames(treatment_joined[3]))
        
        fig1 <- ggplotly(p1) %>% 
          add_annotations(text = paste0("R^2 = ", treatment_r2),
                          x = 0, y = 1, xref = "x", yref = "paper", 
                          xanchor = "left", showarrow = FALSE)
        
        # calculate R2 and plot control replicates
        control_r2 <- round(cor.test(control_joined[[2]],
                                     control_joined[[3]],
                                     method = "pearson")$estimate ^ 2, 2)
        
        p2 <- control_joined %>% 
          ggplot(aes(x = control_joined[[2]], 
                     y = control_joined[[3]],
                     dummy = sgRNA, group = 1)) +
          geom_point() +
          geom_smooth(method = "lm", 
                      formula = "y ~ x",
                      se = FALSE) +
          xlab(colnames(control_joined[2])) +
          ylab(colnames(control_joined[3]))
        
        fig2 <- ggplotly(p2) %>% 
          add_annotations(text = paste0("R^2 = ", control_r2), 
                          x = 0, y = 1, xref = "x", yref = "paper", 
                          xanchor = "left", showarrow = FALSE)
        
        # put the two plots together
        subplot(fig1, fig2,
                titleX = TRUE,
                titleY = TRUE)
      })
      
      ##========================================================================================     
      
      # INDIVIDUAL sgRNAs TAB
      
      # choose dataset for sgRNA tab
      df_sgRNA <- reactive({
        if (input$showdata == "median_norm"){
          df_sgRNA <- get(paste0(selected_screen(), "_median_norm.sgrna_summary.txt"))
        } else {
          df_sgRNA <- get(paste0(selected_screen(), "_control_norm.sgrna_summary.txt"))
        }
      })
      
      # scatter plot of individual sgRNAs for each gene
      output$scatter_sgrna <- renderPlotly({
        
        # suppress error message that shows if no genes selected
        req(input$comp_data_table_rows_selected)
        
        p5 <- df_sgRNA() %>% 
          filter(Gene %in% genes_from_table()) %>% 
          ggplot(aes_string(x = "Gene", y = input$sgrna_y)) +
          geom_jitter(aes(group = Gene, sgRNA = sgrna, fill = Gene),
                      width = 0.1, alpha = 0.6, shape = 21) +
          stat_summary(aes(group = Gene),
                       fun = median, geom = "crossbar",
                       width = 0.5) +
          scale_fill_viridis_d() +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          xlab(NULL) 
        
        ggplotly(p5, tooltip = c("sgRNA", "y"))
        
      })
      
      ##========================================================================================      
      
      # RANKED GENE SUMMARY TAB
      
      ## Two different plots for ranked gene summary: genebar_topn, genebar_selected
      
      # if ranking top N from all genes
      output$genebar_topn <- renderPlotly({
        
        # suppress error message if nothing selected
        req(input$all_screens_rows_selected)
        
        # get collapsed list of NTCs for specific library
        ntc_list <- all_ntc_list %>% 
          filter(library_name == library_name()) %>% 
          pull(NTC) %>% 
          str_c(collapse = "|")
        
        df_top20 <- df_gene() %>%
          select(id, `neg|score`, `neg|rank`,
                 `pos|rank`, `pos|score`,
                 `neg|p-value`, `pos|p-value`,
                 `neg|fdr`, `pos|fdr`) %>%
          pivot_longer(cols = ends_with("rank"), 
                       names_to = "rank_type", values_to = "rank") %>%
          pivot_longer(cols = ends_with("score"), 
                       names_to = "score_type", values_to = "score") %>%
          pivot_longer(cols = ends_with("value"), 
                       names_to = "p_type", values_to = "p_value") %>%
          pivot_longer(cols = ends_with("fdr"), 
                       names_to = "fdr_type", values_to = "fdr") %>%
          arrange(rank) %>%
          filter(rank <= input$topn) %>%
          # only look at neg|score for neg|rank and pos|score for pos_rank
          filter(case_when(rank_type == "neg|rank" ~ score_type == "neg|score",
                           rank_type == "pos|rank" ~ score_type == "pos|score")) %>%
          filter(case_when(rank_type == "neg|rank" ~ p_type == "neg|p-value",
                           rank_type == "pos|rank" ~ p_type == "pos|p-value")) %>% 
          filter(case_when(rank_type == "neg|rank" ~ fdr_type == "neg|fdr",
                           rank_type == "pos|rank" ~ fdr_type == "pos|fdr"))
        
        p3 <- df_top20 %>%
          mutate(is_ntc = case_when(str_detect(ntc_list, id) ~ TRUE,
                                    TRUE ~ FALSE)) %>% 
          ggplot(aes(y = fct_reorder(id, -rank), x = -log10(score),
                     fill = case_when(input$fillby == "ntc" ~ is_ntc,
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
          xlab("-log10 MAGeCK gene score")
        
        gp3 <- ggplotly(p3, tooltip = c("rank", "score", "p", "fdr"))
        
        gp3[["x"]][["layout"]][["annotations"]][[1]][["y"]] <- -0.08
        
        gp3 %>% layout(margin = list(b = 75))
        
      })  
      
      # if plotting only selected genes
      
      output$genebar_selected <- renderPlotly({
        
        # suppress error message if nothing selected
        req(input$all_screens_rows_selected)
        
        ntc_list <- all_ntc_list %>% 
          filter(library_name == library_name()) %>% 
          pull(NTC) %>% 
          str_c(collapse = "|")
        
        df_gene_selected <-  df_gene() %>% 
          filter(id %in% genes_from_table()) %>% 
          select(id, `neg|score`, `neg|rank`, 
                 `pos|rank`, `pos|score`,
                 `neg|p-value`, `pos|p-value`,
                 `neg|fdr`, `pos|fdr`) %>%
          pivot_longer(cols = ends_with("rank"), 
                       names_to = "rank_type", values_to = "rank") %>% 
          pivot_longer(cols = ends_with("score"), 
                       names_to = "score_type", values_to = "score") %>% 
          pivot_longer(cols = ends_with("value"), 
                       names_to = "p_type", values_to = "p_value") %>% 
          pivot_longer(cols = ends_with("fdr"), 
                       names_to = "fdr_type", values_to = "fdr") %>%
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
                 id2 = reorder_within(as.factor(id), -rank, rank_type),
                 is_ntc = case_when(str_detect(ntc_list, id) ~ TRUE,
                                    TRUE ~ FALSE)) %>% 
          ggplot(aes(y = id2, x = -log10(score),
                     fill = case_when(input$fillby == "ntc" ~ is_ntc,
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
          xlab("-log10 MAGeCK gene score") 
        
        gp4 <- ggplotly(p4, tooltip = c("rank", "score", "p", "fdr"))
  
        gp4[["x"]][["layout"]][["annotations"]][[1]][["y"]] <- -0.08
        
        gp4 %>% layout(margin = list(b = 75))
        
      })
      
      ##========================================================================================
      
      # DOT PLOT TAB
      
      output$dotplot <- renderPlotly({
        
        # suppress error message if nothing selected
        req(input$all_screens_rows_selected)
        
        ntc_list <- all_ntc_list %>% 
          filter(library_name == library_name()) %>% 
          pull(NTC) %>% 
          str_c(collapse = "|")
        
        df_gene_is_ntc <- df_gene() %>% 
          mutate(is_ntc = case_when(str_detect(ntc_list, id) ~ TRUE,
                                    TRUE ~ FALSE))
        
        median_neg_ntc <- df_gene_is_ntc %>% 
          filter(is_ntc == TRUE) %>% 
          pull(`neg|score`) %>% 
          median()
        
        median_pos_ntc <- df_gene_is_ntc %>% 
          filter(is_ntc == TRUE) %>% 
          pull(`pos|score`) %>% 
          median()
        
        p7 <- df_gene_is_ntc %>% 
          select(id, is_ntc, `neg|score`, `pos|score`) %>%
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
          geom_point(shape = 21, alpha = 0.6, size = 2) +
          geom_point(data = . %>% filter(is_ntc == TRUE),
                     fill = "white", shape = 21, size = 2) +
          scale_fill_brewer(palette = "Dark2") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.title = element_blank()) +
          xlab("Genes") +
          ylab("Log10 fold change from NTC median")
        
        ggplotly(p7, tooltip = "text")
        
      })
      
      ##========================================================================================
      
      # COMPARE 2 SCREENS TAB
      
      # show all screens again, this time with multiple selection
      output$all_screens_2 <- renderDataTable({
        metadata %>% 
          datatable(rownames = FALSE,
                    filter = "top",
                    extensions = c("Buttons", "ColReorder"),
                    options = list(search = list(regex = TRUE),
                                   dom = "Bfrtip",
                                   buttons = I("colvis"),
                                   colReorder = list(realtime = FALSE)),
                    escape = FALSE) 
      })
      
      # get screen names from selected datatable rows
      screen_choices <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_2_rows_selected) %>%  
          pull(name)
      })
      
      # get library names from selected datatable rows
      library_choices <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_2_rows_selected) %>% 
          pull(LibraryName)
      })
      
      # print selected screen names in sidebar
      output$selected_2_screens <- renderText({
        HTML(paste(screen_choices()[1], screen_choices()[2], sep = "<br>"))
      })
      
      # choose output file type to use
      output_file_type <- reactive({
        if (input$norm_type == "median_norm"){
          output_file_type <- "median_norm.gene_summary.txt"
        } else {
          output_file_type <- "control_norm.gene_summary.txt"
        }
      })
      
      # assemble plot data
      compare_scores_ranks <- reactive({
        
        # get NTC list for each screen's library
        ntc_list_screen1 <- all_ntc_list %>% 
          filter(library_name == library_choices()[1]) %>% 
          pull(NTC) %>% 
          str_c(collapse = "|")
        
        ntc_list_screen2 <- all_ntc_list %>% 
          filter(library_name == library_choices()[2]) %>% 
          pull(NTC) %>% 
          str_c(collapse = "|")
        
        # set output df names
        # and add is_ntc column
        screen1_output_df <- get(paste0(screen_choices()[1], "_", output_file_type())) %>% 
          mutate(is_ntc = case_when(str_detect(ntc_list_screen1, id) ~ TRUE,
                                    TRUE ~ FALSE))
        
        screen2_output_df <- get(paste0(screen_choices()[2], "_", output_file_type())) %>% 
          mutate(is_ntc = case_when(str_detect(ntc_list_screen2, id) ~ TRUE,
                                    TRUE ~ FALSE))
        
        # if libraries from both screens are different, remove all NTCs
        # if libraries are same, include NTCs
        if (library_choices()[1] != library_choices()[2]){
          
          screen1_output_df <- screen1_output_df %>% 
            filter(is_ntc == FALSE)
          
          screen2_output_df <- screen2_output_df %>% 
            filter(is_ntc == FALSE)
          
        }
        
        # only comparing genes in common between the 2 screens
        compare_scores <- screen1_output_df %>% 
          select(id, `neg|score`, `pos|score`) %>% 
          full_join(select(screen2_output_df, 
                           id, `neg|score`, `pos|score`),
                    by = "id",
                    suffix = c("_screen1", "_screen2")) %>% 
          pivot_longer(cols = contains("score"), 
                       names_to = "score_type", values_to = "score") %>% 
          separate(score_type, into = c("score_type", "screen"), sep = "_")
        
        compare_ranks <- screen1_output_df %>% 
          select(id, `neg|rank`, `pos|rank`) %>% 
          full_join(select(screen2_output_df, 
                           id, `neg|rank`, `pos|rank`),
                    by = "id",
                    suffix = c("_screen1", "_screen2")) %>% 
          pivot_longer(cols = contains("rank"), 
                       names_to = "rank_type", values_to = "rank") %>% 
          separate(rank_type, into = c("rank_type", "screen"), sep = "_")
        
        # join ranks & scores
        # add key for plotly_click
        compare_scores_ranks <- full_join(compare_scores, compare_ranks, 
                                          by = c("id", "screen")) %>% 
          filter((str_detect(score_type, "neg") & str_detect(rank_type, "neg")) |
                   (str_detect(score_type, "pos") & str_detect(rank_type, "pos"))) %>% 
          mutate(score_type = str_replace(score_type, "\\|.*$", "")) %>% 
          select(-rank_type) %>% 
          pivot_wider(names_from = screen, values_from = c(score, rank)) %>% 
          mutate(key = row.names(.))
        
      })
      
      # plot 2 screens
      output$compare_plot <- renderPlotly({
        
        # suppress error message if nothing selected
        req(input$all_screens_2_rows_selected)
        
        # make this df non-reactive, as reactive seems to generate an error with slice_min 
        compare_scores_ranks <- compare_scores_ranks()
        
        
        p8 <- compare_scores_ranks %>% 
          ggplot(aes(x = -log10(score_screen1), y = -log10(score_screen2), 
                     id = id, key = key, 
                     rank_screen1 = rank_screen1, rank_screen2 = rank_screen2)) +
          geom_point(shape = 21, alpha = 0.5, size = 2) +
          geom_abline(slope = 1, linetype = "dashed") +
          facet_wrap(~score_type, scales = "free") +
          xlab(paste0("\n", screen_choices()[1], "(screen1)")) +
          ylab(paste0(screen_choices()[2], "\n(screen2)")) 
        
        ggplotly(p8, source = "compare2_plot",
                 tooltip = c("id", "x", "y", "rank_screen1", "rank_screen2")) %>% 
          layout(margin = list(l = 75, b = 75))
        
        gp8 <- ggplotly(p8, source = "compare2_plot",
                       tooltip = c("id", "x", "y", "rank_screen1", "rank_screen2")) 
        
        # move the axis titles so they're not squished against the tick labels
        gp8[["x"]][["layout"]][["annotations"]][[2]][["x"]] <- -0.05
        gp8[["x"]][["layout"]][["annotations"]][[1]][["y"]] <- -0.08
        
        gp8 %>% layout(margin = list(l = 75, b = 75))
        
      })
      
      # if point on plot is clicked, show data for that gene across all screens
      output$gene_all_screens_table <- renderDataTable({
        
        #suppress error message
        req(event_data(source = "compare2_plot",
                       "plotly_click"))
        
        # get gene name for clicked point
        click_data_key <- event_data(source = "compare2_plot",
                                     "plotly_click") %>% 
          pull(key)
        
        click_gene <- compare_scores_ranks() %>% 
          filter(key == click_data_key) %>% 
          pull(id)
        
        # get all-screens joined dataset
        # and filter by gene name
        all_screens_output_df <- get(paste0(output_file_type(), "_joined"))
        
        all_screens_output_df %>% 
          filter(id == click_gene) %>% 
          datatable(rownames = FALSE,
                    filter = "top",
                    extensions = c("Buttons", "ColReorder"),
                    options = list(search = list(regex = TRUE),
                                   dom = "Bfrtip",
                                   buttons = I("colvis"),
                                   colReorder = list(realtime = FALSE)),
                    escape = FALSE)
        
      })
      
    }
    
  })
  
  
}

shinyApp(ui = ui, server = server)
