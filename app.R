
library(shiny)
library(synapser)
library(synapserutils)
library(tidyverse)
library(tidytext)
library(DT)
library(plotly)
library(waiter)

ui <- fluidPage(
  
  # Synapse stuff from template
  tags$head(
    singleton(
      includeScript("www/readCookie.js")
    )
  ),
  
  titlePanel("HIV-CRISPR data viz demo (CUL3 library only)"),
  
  uiOutput("title"),
  
  tabsetPanel(
    tabPanel("Metadata",
             sidebarLayout(
               sidebarPanel(
                 helpText("This table shows metadata for all screens. 
                        Select a screen to use in the next 5 tabs.",
                        br(), br(),
                        "Drag columns to reorder, 
                        or click 'Column visibility' to show/hide columns. 
                        Click to select multiple rows for plotting. 
                        Selected genes will be listed below. 
                        Search box accepts regular expressions (ie. use '|' for OR).",
                        br(), br()),
                 p("Selected screen: ",
                   textOutput("show_selected_screen", inline = TRUE))
               ),
               mainPanel(
                 br(),
                 dataTableOutput("all_screens"),
                 hr(),
                 dataTableOutput("single_metadata")
               )
             )
             
    ),
    tabPanel("Output data",
             sidebarLayout(
               sidebarPanel(
                 helpText("This table shows output data for the selected screen. 
                        Select genes to plot in the next 3 tabs.",
                        br(), br(),
                        "Drag columns to reorder, 
                          or click 'Column visibility' to show/hide columns. 
                          Click to select multiple rows for plotting. 
                          Selected genes will be listed below. 
                          Search box accepts regular expressions (ie. use '|' for OR).",
                        br(), br()),
                 radioButtons("showdata", "Select output files to show:",
                              choices = list("median_norm.gene_summary" = "median_norm",
                                             "control_norm.gene_summary" = "control_norm"),
                              selected = "median_norm"),
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
                 helpText("This plot only works with median_norm.gene_summary.txt")
               ),
               mainPanel(
                 plotlyOutput("dotplot")
               )
             )
             
    ),
    tabPanel("Compare 2 Screens",
             sidebarLayout(
               sidebarPanel(
                 p("Selected screens: ",
                   br(),
                   htmlOutput("selected_2_screens", inline = TRUE)),
                 radioButtons("norm_type", "Select output files to use:",
                              choices = list("median_norm.gene_summary" = "median_norm",
                                             "control_norm.gene_summary" = "control_norm"),
                              selected = "median_norm"),
               ),
               mainPanel(
                 br(),
                 dataTableOutput("all_screens_2"),
                 hr(),
                 plotlyOutput("compare_plot"),
                 br(),
                 dataTableOutput("gene_all_screens_table")
               )
             )
             
    )
  ),
  
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
      ### login and update session; otherwise, notify to login to Synapse first
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
              "There was an error with the login process. Please refresh your Synapse session by logging out of and back in to",
              a("Synapse", href = "https://www.synapse.org/", target = "_blank"),
              ", then refresh this page."
            )
          )
        )
        
      })
      
      # Any shiny app functionality that uses synapse should be within the
      # input$cookie observer (so basically everything)
      output$title <- renderUI({
        titlePanel(sprintf("Welcome, %s", synGetUserProfile()$userName))
        
      })
      
      # Load data from Synapse
      source("get_synapse_data.R")
      
      # Show metadata for all screens that have a configId
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
      
      # get synID for chosen screen
      sample_sheet_id <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_rows_selected) %>%  
          pull(configId)
      })
      
      # get selected screen id and then print (separate bc need to use selected_screen again later)
      selected_screen <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_rows_selected) %>%  
          pull(name)
      })
      
      # print selected screen in sidebar
      output$show_selected_screen <- renderText({
        selected_screen()
      })
      
      # pull library name for eventual use in identifying NTCs
      library_name <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_rows_selected) %>%  
          pull(LibraryName)        
      })
      
      # single screen metadata (from treatment replicate 1)
      output$single_metadata <- renderDataTable({
        
        # get data for specified screen
        get_counts(sample_sheet_id())
        
        metadata %>% 
          filter(configId == sample_sheet_id()) %>%
          select(-starts_with("ROW_")) %>% 
          pivot_longer(cols = everything(),
                       names_to = "field", values_to = "value",
                       values_transform = list(value = as.character)) %>% 
          datatable(caption = paste0("Metadata for screen ", screen_name, ", pulled from treatment replicate 1"),
                    rownames = FALSE)
      })
      
      # choose dataset for Data tab
      df_gene <- reactive({
        if (input$showdata == "median_norm"){
          df_gene <- get(paste0(selected_screen(), "_median_norm.gene_summary.txt"))
        } else {
          df_gene <- get(paste0(selected_screen(), "_control_norm.gene_summary.txt"))
        }
      })
      
      # show data for single screen
      output$comp_data_table <- renderDataTable({
        
        # bring in GeneCards links
        # don't replace gene id with link - there are a few differences between `MyList` and `Gene Symbol`
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
      
      # print list of selected genes in sidebar
      # do it multiple times bc you can't reuse same output
      output$selection_info1 <- 
        output$selection_info2 <- 
        output$selection_info3 <-
        output$selection_info4 <-
        renderText({
          genes_from_table()
        },
        sep = ", ")
      
      # QC scatter plot of count files with R2 values
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
                titleY = TRUE) %>% 
          layout(title = "Counts in treatment (left) and control (right) replicates",
                 margin = list(t = 50))
      })
      
      
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
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          xlab(NULL) +
          ggtitle("Individual sgRNA scores for selected genes")
        
        ggplotly(p5, tooltip = c("sgRNA", "y"))
        
      })
      
      
      ## Two different plots for ranked gene summary: genebar_topn, genebar_selected
      
      # if ranking top N from all genes
      output$genebar_topn <- renderPlotly({
        
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
          pivot_longer(cols = ends_with("rank"), names_to = "rank_type", values_to = "rank") %>%
          pivot_longer(cols = ends_with("score"), names_to = "score_type", values_to = "score") %>%
          pivot_longer(cols = ends_with("value"), names_to = "p_type", values_to = "p_value") %>%
          pivot_longer(cols = ends_with("fdr"), names_to = "fdr_type", values_to = "fdr") %>%
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
          xlab("-log10 MAGeCK gene score") +
          ggtitle(paste0("Top ", input$topn, " highest-ranked genes in ", input$showdata, ".gene_summary"))
        
        ggplotly(p3, tooltip = c("rank", "score", "p", "fdr"))
        
      })  
      
      output$genebar_selected <- renderPlotly({
        
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
          xlab("-log10 MAGeCK gene score") +
          ggtitle(paste0("Selected genes from ", input$showdata, ".gene_summary, by pos or neg rank"))
        
        ggplotly(p4, tooltip = c("rank", "score", "p", "fdr"))
        
      })
      
      output$dotplot <- renderPlotly({
        
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
          geom_point(shape = 21, alpha = 0.8, size = 2) +
          geom_point(data = . %>% filter(is_ntc == TRUE),
                     shape = 21, fill = "white", size = 2, 
                     show.legend = FALSE) +
          scale_fill_brewer(palette = "Dark2") +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) +
          xlab("Genes") +
          ylab("Log10 fold change from NTC median")
        
        ggplotly(p7, tooltip = "text")
        
      })
      
      # show all screens again, this time with multiple selection
      output$all_screens_2 <- renderDataTable({
        metadata %>% 
          datatable(rownames = FALSE) 
      })
      
      # get screen names from selected datatable rows
      screen_choices <- reactive({
        metadata %>% 
          filter(row_number() %in% input$all_screens_2_rows_selected) %>%  
          pull(name)
      })
      
      # print selected screen names
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
      
      compare_scores_ranks <- reactive({
        
        # get screen names for each selected screen
        # screen1 <- screen_choices()[1]
        # screen2 <- screen_choices()[2]
        
        # set output df names
        screen1_output_df <- get(paste0(screen_choices()[1], "_", output_file_type()))
        screen2_output_df <- get(paste0(screen_choices()[2], "_", output_file_type()))
        
        # only comparing genes in common between the 2 screens
        compare_scores <- screen1_output_df %>% 
          select(id, `neg|score`, `pos|score`) %>% 
          full_join(select(screen2_output_df, 
                           id, `neg|score`, `pos|score`),
                    by = "id",
                    suffix = c("_screen1", "_screen2")) %>% 
          pivot_longer(cols = contains("score"), names_to = "score_type", values_to = "score") %>% 
          separate(score_type, into = c("score_type", "screen"), sep = "_")
        
        compare_ranks <- screen1_output_df %>% 
          select(id, `neg|rank`, `pos|rank`) %>% 
          full_join(select(screen2_output_df, 
                           id, `neg|rank`, `pos|rank`),
                    by = "id",
                    suffix = c("_screen1", "_screen2")) %>% 
          pivot_longer(cols = contains("rank"), names_to = "rank_type", values_to = "rank") %>% 
          separate(rank_type, into = c("rank_type", "screen"), sep = "_")
        
        # join ranks & scores
        # add key for plotly_click
        compare_scores_ranks <- full_join(compare_scores, compare_ranks, by = c("id", "screen")) %>% 
          filter((str_detect(score_type, "neg") & str_detect(rank_type, "neg")) |
                   (str_detect(score_type, "pos") & str_detect(rank_type, "pos"))) %>% 
          mutate(score_type = str_replace(score_type, "\\|.*$", "")) %>% 
          select(-rank_type) %>% 
          pivot_wider(names_from = screen, values_from = c(score, rank)) %>% 
          mutate(key = row.names(.))
        
        # set up key for plotly_click
        #  compare_scores_ranks$key <- row.names(compare_scores_ranks)
        
      })
      
      output$compare_plot <- renderPlotly({
        
        # make this df non-reactive since reactive seems to generate an error with slice_min 
        compare_scores_ranks <- compare_scores_ranks()
        
        p8 <- compare_scores_ranks %>% 
          ggplot(aes(x = -log10(score_screen1), y = -log10(score_screen2), 
                     id = id, key = key, rank_screen1 = rank_screen1, rank_screen2 = rank_screen2)) +
          geom_point(shape = 21, alpha = 0.5) +
          geom_abline(slope = 1) +
          facet_wrap(~score_type, scales = "free") +
          geom_point(data = compare_scores_ranks %>% slice_min(rank_screen1, n = 20),
                     shape = 21, size = 3, fill = "turquoise4") +
          geom_point(data = compare_scores_ranks %>% slice_min(rank_screen2, n = 20),
                     shape = 21, size = 3, fill = "turquoise4") +
          xlab(paste0(screen_choices()[1], "\n(screen1)")) +
          ylab(paste0(screen_choices()[2], "\n(screen2)")) +
          ggtitle("Comparison of -log10(score) for 2 screens, with top 20 hits for each highlighted")
        
        ggplotly(p8, source = "compare2_plot",
                 tooltip = c("id", "x", "y", "rank_screen1", "rank_screen2")) %>% 
          layout(margin = list(l = 150))
        
      })
      
      output$gene_all_screens_table <- renderDataTable({
        
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
          datatable()
        
      })
      
    }
    
  })
  
  
}

shinyApp(ui = ui, server = server)
