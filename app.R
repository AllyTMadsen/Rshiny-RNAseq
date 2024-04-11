#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(BiocManager)
library(SummarizedExperiment)
library(plotly)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(fields)
library(fgsea)
library(ggbeeswarm)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(

  theme = bs_theme(bootswatch = "journal"),  
  titlePanel("BF591 Final Project: Allison Madsen"),
  
  # Define tabs with different sidebar and main panel for each
  tabsetPanel(
    
    # Tab 1
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 h4("Upload Raw Counts Data to Start Analysis: "),
                 fileInput("file1", label = NULL, placeholder = ".csv, .tsv, .txt", accept = c('.tsv', '.csv', '.txt')),
                 p("If over 75 unique sample IDs are detected, use the slider to select a number of samples to display:"),
                 sliderInput("sample_number_slider", "Samples to display:", min = 0, max = 75, value = 12),
               ),
               mainPanel(
                 # Content for tab 1 goes here
                 tabsetPanel(
                   
                   # Sub-tab 1
                   tabPanel("Data Summary",
                            # Content for Sub-tab 1.1 goes here
                            h4("Table Including the Type and Summary of Values in each Column: "),
                            tableOutput("tab1_summary_table")
                   ),
                   
                   # Sub-tab 2
                   tabPanel("Sample Information",
                            # Content for Sub-tab 1.2 goes here
                            h4("Sortable Table of Sample Information for each Gene"),
                            DTOutput("tab1_sample_info_table")
                   ),
                   
                   # Sub-tab 3
                   tabPanel("Plot of Continuous Variables",
                            # Content for Sub-tab 1.3 goes here
                            h4("Histogram for Total Number of Raw Counts across Different Samples"),
                            plotOutput('tab1_continuous_plots')
                   )
                 )
               )
             )
    ),
    
    # Tab 2
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 h3("Use the Button Below to Normalize Provided Counts Data: "),
                 actionButton("run_DE_norm", label = "Run DESeq Normalization!"),
                 p("Use the sliders to filter the data based on Variance and Number of Non-zero Genes: "),
                 sliderInput("var_slider", "% Variance to include: ", min = 0, max = 100, value = 25),
                 sliderInput("zero_slider", "Number of non-zero genes permitted per sample: ", min = 0, max = 7, value = 5),
                 p("Use the Button Below to Update the Plots after Adjusting the Slider Values: "),
                 actionButton("update_counts_plots", "Update Plots!")
               ),
               mainPanel(
                 # Content for tab 2 goes here
                 tabsetPanel(
                   # Sub-tab 1
                   tabPanel("Filter Summary",
                            h3("Filter Summary of Normalized Counts data: "),
                            tableOutput("Filter_summary_table")
                   ),
                   # Sub-tab 2
                   tabPanel("Diagnostic Scatter Plots",
                            h4("Plots for Median Count vs Variance and Number of Zeros"),
                            p("Samples that pass the filter are displayed in dark blue color.  Samples that have been filtered out are lighter"),
                            plotOutput("scatter_plot1"),
                            plotOutput("scatter_plot2")
                            
                   ),
                   # Sub-tab 3
                   tabPanel("Heatmap",
                            h4("Clustered Heatmap of Filtered Gene Samples:"),
                            plotOutput("clustered_heatmap")
                   ),
                   # Sub-tab 4
                   tabPanel("Interactive Scatter Plot of PCA Projections",
                            fluidRow(
                              column(3,
                                     h4("Select PCs to plot with buttons below: "),
                                     # Add sidebar content here
                                     radioButtons("pc_selection1", "Select Principal Components:", choices = 1:5, selected = 1),
                                     radioButtons("pc_selection2", "Select Principal Components:", choices = 1:5, selected = 2)
                              ),
                              column(9,
                                     h3("PCA plot for selected PCs:"),
                                     plotOutput("PCA_tab2")
                              )
                            )
                   )
                 )
               )
             )
    ),
    
    # Tab 3
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 h3("Upload Meta Data File Below: "),
                 fileInput("file2", label = NULL, placeholder = ".csv, ,tsv, .txt", accept = c('.tsv', '.csv', '.txt')),
                 h4("Press Button Below to Run DESeq2: "),
                 actionButton("run_DE_dds", label = "Run DESeq DE analysis!"),
                 #sliderInput("dds_gene_num", "Number of rows to display: ", min = 0, max = 32613, value = 50),
               ),
               mainPanel(
                 # Content for tab 3 goes here
                 tabsetPanel(
                   # Sub-tab 1
                   tabPanel("Sortable Table",
                            h4("Table of DESeq2 Results wtih Sortable Columns"),
                            DTOutput("DESEQ_res_table")
                   ),
                   # Sub-tab 2-- removed, code still exists
                   # tabPanel("DE Assignment Volcano Plot",
                   #          h4("Volcano Plot of Differential Expression Results"),
                   #          sliderInput("padj_threshold", "padj threshold to include: ", min = 0, max = 1.0, value = 0.1),
                   #          plotOutput("Volcano_plot")
                   # ),
                   tabPanel("R Shiny Assignment Volcano Plot",
                     h4("A volcano plot can be generated with 'log2 fold-change' on the x-axis and 'p-adjusted' on the y-axis."),
                     fluidRow(
                       column(4, 
                              radioButtons("button_x", "Choose the column for the x-axis",
                                           choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                              radioButtons("button_y", "Choose the column for the y-axis",
                                           choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                              colourInput("color1", "Base Point Color", "maroon"),
                              colourInput("color2", "Highlight Point Color", "orange"),
                              sliderInput(inputId = "padj_magnitude", min = -300, max = 0, value = -150, label = "Select the magnitude of the p adjusted coloring:"),
                              actionButton("go", label = "Plot", class = "btn-primary icon-button", icon = icon("fire"), width = "100%")
                              ),
                       column(8,
                              plotOutput("rshiny_plot")
                              )
                     )
                   )
                 )
               )
             )
    ),
    #tab 4
    tabPanel("Individual Gene Expression",
             sidebarLayout(
               sidebarPanel(
                 h3("Choose a Catagorical Field and Gene ID to Visualize Individual Expression Data."),
                 radioButtons("meta_choice", "Please select metadata information: ", choices = c("replicate", "sample", "timepoint"), selected = "timepoint"),
                 radioButtons("plot_choice", "Please select desired plot type: ", choices = c("box plot", "violin plot", "beeswarm plot", "bar plot"), selected = "bar plot"),
                 textInput("gene_search", "Search for Gene ID:"),
                 actionButton("go_choice", label = "click to generate plot!")
               ),
               mainPanel(
                 h4("Plots of Normalized Counts Data for Specified Gene"),
                 plotOutput("choice_plot"),
               )
             ))
  )
  )


# Define server logic 
server <- function(input, output, session) {
  
  norm_button_clicked <- reactiveVal(FALSE)
  norm_counts <- reactiveVal(NULL)
  result_list <- reactiveValues(results_df = NULL, dds = NULL)
  filtered_results <- reactiveVal(NULL)
  fgsea_results <- reactiveVal(NULL)
  plot_list <- reactiveValues()
  last_plot <- reactiveVal(NULL)
  go_button_clicked <- reactiveVal(FALSE)
  
  
  
  #create raw data reactive val-- this time filter out genes with 0 variance
  raw_data <- reactive({
    req(input$file1)
    data1 <- read.table(input$file1$datapath, header = TRUE, sep = "\t")
    
    nums <- data1[, -1, drop = FALSE]
    row_var <- apply(nums, 1, var, na.rm = TRUE)
    nonzero_var_rows <- which(row_var != 0)
    return(data1[nonzero_var_rows, ])
  })
  
  #output summary table tab 1
  output$tab1_summary_table <- renderTable({
    req(input$file1)
    
    summary_data <- data.frame(
      Column_Names = names(raw_data()),
      Type = sapply(raw_data(), function(x) class(x)),
      Mean_stdev_Distinct_Values = sapply(raw_data(), function(x) {
        if (is.numeric(x)) {
          mean_val <- ifelse(is.integer(x), as.integer(mean(x, na.rm = TRUE)), mean(x, na.rm = TRUE))
          paste(mean_val, "+/- ", round(sd(x, na.rm = TRUE), 2))
        } else {
          unique_vals <- unique(x)
          if (length(unique_vals) > 100) {
            unique_vals <- unique_vals[1:input$sample_number_slider]
          }
          paste(unique_vals, collapse = ", ")
        }
      }),
      stringsAsFactors = FALSE
    )
    
    return(summary_data)
  })
  
  #tab 1 sample info table
  output$tab1_sample_info_table <- renderDT({
    req(input$file1)
    return(raw_data())
  })
  
  
  #tab 1 subtab 3 counst histogram data
  output$tab1_continuous_plots <- renderPlot({
    req(input$file1, raw_data())
    
    # Select only numeric columns for violin plots
    numeric_cols <- Filter(is.numeric, raw_data())
    
    data_long <- tidyr::gather(numeric_cols, key = "sample_name", value = "count", everything())
    
    ggplot(data_long, aes(x = sample_name, y = count, fill = sample_name)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
      labs(title = "Histogram of Counts by Sample",
           x = "Sample Name",
           y = "Count",
           fill = "Sample Name")
  })
    
  #run counts normalization
  observeEvent(input$run_DE_norm, {
    # Run DESeq normalization when the button is clicked
    req(input$file1)
    
    count_mat <- as.matrix(raw_data()[-1])
    row.names(count_mat) <- raw_data()$gene
    
    dds <- DESeqDataSetFromMatrix(
      countData = count_mat,
      colData = data.frame(sample_name = colnames(raw_data()[-1])),
      design = ~1
    )
    
    dds <- estimateSizeFactors(dds)
    
    # Store the normalized counts for later use
    norm_counts(data.frame(fpm(dds, robust = FALSE)) %>%    #this was a tibble using as_tibble
                  mutate(gene = raw_data()$gene) %>%
                  relocate(gene))
  })
  
  
  output$Normalized_counts_table <- renderTable({
    req(input$file1, norm_counts())
    head(norm_counts(), 20)
  })
  
  
  # Reactive expression for filtered genes
  filtered_genes <- reactive({
    req(norm_counts())
    
    # Filtering based on var_slider
    var_filter <- apply(as.matrix(norm_counts()[-1]), 1, sd) >= (input$var_slider / 100)
    
    # Filtering based on zero_slider
    zero_filter <- apply(as.matrix(norm_counts()[-1] > 0), 2, sum) >= input$zero_slider
    
    # Apply the filters
    filtered_genes <- which(zero_filter & var_filter)
    
    return(filtered_genes)
  })
  
  # Reactive expression for filtered counts
  filtered_counts <- reactive({
    req(filtered_genes(), norm_counts())
    return(norm_counts()[filtered_genes(), ])
  })
  
  
  #FUNCTION
  
  #function for summary table tab 2 subtab1
  generateFilterSummaryTable <- function(norm_counts, var_slider, zero_slider) {
    # Filtering based on var_slider
    var_threshold <- apply(as.matrix(norm_counts()[-1]), 1, var) >= (input$var_slider / 100)
    passed_var_filter <- sum(var_threshold)
    
    # Filtering based on zero_slider
    passed_zero_filter <- sum(apply(as.matrix(t(norm_counts()[-1])), 2, sum) >= input$zero_slider)
  
    # Calculate summary values
    total_genes <- nrow(norm_counts())
    
    filter_summary_data <- data.frame(
      Metric = c("Number of Samples", "Total Number of Genes", "Number of Genes Passing Var Filter", "Percentage of Genes Passing Var Filter",
                 "Number of Genes Passing Zero Filter", "Percentage of Genes Passing Zero Filter"),
      Value = c(ncol(norm_counts()), total_genes, passed_var_filter, passed_var_filter / total_genes * 100,
                passed_zero_filter, passed_zero_filter / total_genes * 100)
    )
    
    return(filter_summary_data)
  }
  
  # Render the filter summary table
  output$Filter_summary_table <- renderTable({
    req(norm_counts())
    generateFilterSummaryTable(norm_counts(), input$var_slider, input$zero_slider)
  })
  
  
  # Create diagnostic scatter plots
  
  #scatter plot for mean vs variance
  output$scatter_plot1 <- renderPlot({
    req(norm_counts())
    
    median_exp <- apply(norm_counts()[-1], 1, median)
    gene_var <- apply(norm_counts()[-1], 1, var)
    
    var_threshold <- apply(as.matrix(norm_counts()[-1]), 1, var) >= (input$var_slider / 100)
    zero_threshold <- apply(as.matrix(t(norm_counts()[-1])), 2, sum) >= input$zero_slider
    
    plot_data <- data.frame(Median_expression = median_exp,
                            Variance = gene_var,
                            PassedVarFilter = var_threshold,
                            PassedZeroFilter = zero_threshold)
    
    plot_data$Color <- ifelse(plot_data$PassedVarFilter & plot_data$PassedZeroFilter, "darkblue", "lightblue")
    
    plot <- ggplot(plot_data, aes(x = Variance, y = Median_expression, color = Color)) + 
      geom_point(aes(alpha = 0.5), size = 1) +
      xlab("Variance (log10 scale)") + 
      ylab("Median Count") +
      ggtitle("Median Count/Variance Plot") +
      theme(legend.position = "none") +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(values = c("darkblue", "lightblue"))
    
    return(plot)
  })
  
  
  
  # Scatter plot for mean_exp vs total number of zeros
  output$scatter_plot2 <- renderPlot({
    req(norm_counts(), input$var_slider, input$zero_slider)
    
    median_exp <- apply(norm_counts()[-1], 1, median)
    num_zeros <- apply(as.matrix(norm_counts()[-1] == 0), 1, sum)
    
    var_threshold <- apply(as.matrix(norm_counts()[-1]), 1, var) >= (input$var_slider / 100)
    zero_threshold <- apply(as.matrix(t(norm_counts()[-1])), 2, sum) >= input$zero_slider
    
    
    plot_data <- data.frame(Median_expression = median_exp,
                            Variance = num_zeros,
                            PassedVarFilter = var_threshold,
                            PassedZeroFilter = zero_threshold)
    
    plot_data$Color <- ifelse(plot_data$PassedVarFilter & plot_data$PassedZeroFilter, "darkblue", "lightblue")
    
    plot <- ggplot(plot_data, aes(x = Variance, y = Median_expression, color = Color)) + 
      geom_point(aes(alpha = 0.5), size = 1) +
      xlab("Zero Count per Sample") + 
      ylab("Median Count") +
      ggtitle("Median Count vs Number of Zeros Plot") +
      theme(legend.position = "none") +
      scale_x_log10() + scale_y_log10() +
      scale_color_manual(values = c("darkblue", "lightblue"))
    
    
    return(plot)
  })
  
  
  #HEATMAP FUNCTION
  output$clustered_heatmap <- renderPlot({
    req(norm_counts(), input$var_slider, input$zero_slider)
    heatmap_data <- log2(as.matrix(filtered_counts()[-1]) + 1)
    
    # Create the heatmap
    heatmap_obj <- heatmap(heatmap_data, Rowv = NA, scale = "row", labRow = filtered_counts()$gene)
    
    # Add a color bar to the legend
    legend_color_bar <- seq(min(heatmap_data), max(heatmap_data), length.out = 100)
    image.plot(x = rep(1, 100), y = legend_color_bar, z = matrix(legend_color_bar, nrow = 1), 
               col = heat.colors(100), axes = FALSE, legend.only = TRUE, legend.width = 1.5)
    
    # Combine the heatmap and color bar
    layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1.5), heights = c(1, 4))
    return(heatmap_obj)
  })
  
  #FUCNTION TO HELP CREATE se OBJECT
  timepoint_from_sample <- function(x) {
    find <- str_locate(x, "_") 
    timepoint <- str_sub(x, 1, find -1)
    timepoint <- unique(timepoint)
    return(factor(timepoint))
  }
  
  #HELPER FUNCTION TO CREATE META DATA
  sample_replicate <- function(x) {
    x <- sub("_", "", x)  
    x <- str_sub(x, 4, -1)
    return(x)
  }
  
  #PCA PLOT WITH 2 PCs for tab 2 
  output$PCA_tab2 <- renderPlot({
    req(raw_data(), input$var_slider, input$zero_slider)
   
    timepoints <- sapply(colnames(raw_data()[-1]), timepoint_from_sample)
    
    replicates <- sapply(colnames(raw_data()[-1]), sample_replicate)
    #print(replicates)
    meta <- tibble(sample = colnames(raw_data()[-1]), timepoint = timepoints, replicate = replicates)
  
    dds <- DESeqDataSetFromMatrix(
      countData = as.matrix(raw_data()[-1]),
      colData = as.data.frame(meta),
      design = ~1
    )
    
    pca_res <- prcomp(t(assay(dds)))
    pca_data <- as.data.frame(pca_res$x)
    pca_data$timepoint <- meta$timepoint
    
    pc_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
    pc_var_labels <- paste0("PC", 1:ncol(pca_res$x), ": ", round(pc_var, 0), "% Variance")
    
    selected_pc_x <- paste0("PC", input$pc_selection1)
    selected_pc_y <- paste0("PC", input$pc_selection2)
    
    plot <- ggplot(pca_data, aes_string(x = selected_pc_x, y = selected_pc_y, color = "timepoint")) +
      geom_point() +
      labs(title = paste("PCA Plot: Principal Components", selected_pc_x, "vs", selected_pc_y),
           x = paste0(selected_pc_x, ": ", round(pc_var[as.numeric(input$pc_selection1)], 0), "% Variance"),
           y = paste0(selected_pc_y, ": ", round(pc_var[as.numeric(input$pc_selection2)], 0), "% Variance"))
    
    return(plot)
  })
  
  #UPDATING TAB 2 PLOTS AFTER SLIDER UPDATES
  observeEvent(input$update_counts_plots, {
    updateCountsPlots(input$var_slider, input$zero_slider)
  })
  
  #UPdate plots-- all 4 output functions from tab 2 contained below:
  updateCountsPlots <- function(var_slider_val, zero_slider_val) {
    
    output$scatter_plot1 <- renderPlot({
      req(norm_counts())
      
      median_exp <- apply(norm_counts()[-1], 1, median)
      gene_var <- apply(norm_counts()[-1], 1, var)
      
      var_threshold <- apply(as.matrix(norm_counts()[-1]), 1, var) >= (var_slider_val / 100)
      zero_threshold <- apply(as.matrix(t(norm_counts()[-1])), 2, sum) >= zero_slider_val
      
      plot_data <- data.frame(Median_expression = median_exp,
                              Variance = gene_var,
                              PassedVarFilter = var_threshold,
                              PassedZeroFilter = zero_threshold)
      
      plot_data$Color <- ifelse(plot_data$PassedVarFilter & plot_data$PassedZeroFilter, "darkblue", "lightblue")
      
      plot <- ggplot(plot_data, aes(x = Variance, y = Median_expression, color = Color)) + 
        geom_point(aes(alpha = 0.5), size = 1) +
        xlab("Variance (log10 scale)") + 
        ylab("Median Count") +
        ggtitle("Median Count/Variance Plot") +
        theme(legend.position = "none") +
        scale_x_log10() + scale_y_log10() +
        scale_color_manual(values = c("darkblue", "lightblue"))
      
      return(plot)
    })
    
    output$scatter_plot2 <- renderPlot({
      req(norm_counts())
      
      median_exp <- apply(norm_counts()[-1], 1, median)
      num_zeros <- apply(as.matrix(norm_counts()[-1] == 0), 1, sum)
      
      var_threshold <- apply(as.matrix(norm_counts()[-1]), 1, var) >= (var_slider_val / 100)
      zero_threshold <- apply(as.matrix(t(norm_counts()[-1])), 2, sum) >= zero_slider_val
      
      
      plot_data <- data.frame(Median_expression = median_exp,
                              Variance = num_zeros,
                              PassedVarFilter = var_threshold,
                              PassedZeroFilter = zero_threshold)
      
      plot_data$Color <- ifelse(plot_data$PassedVarFilter & plot_data$PassedZeroFilter, "darkblue", "lightblue")
      
      plot <- ggplot(plot_data, aes(x = Variance, y = Median_expression, color = Color)) + 
        geom_point(aes(alpha = 0.5), size = 1) +
        xlab("Zero Count per Sample") + 
        ylab("Median Count") +
        ggtitle("Median Count vs Number of Zeros Plot") +
        theme(legend.position = "none") +
        scale_x_log10() + scale_y_log10() +
        scale_color_manual(values = c("darkblue", "lightblue"))
      
      
      return(plot)
    })
    
    output$clustered_heatmap <- renderPlot({
      heatmap_data <- log2(as.matrix(filtered_counts()[-1]) + 1)
      #print(filtered_counts())
      
      # Create the heatmap
      heatmap_obj <- heatmap(heatmap_data, Rowv = NA, scale = "row", labRow = filtered_counts()$gene)
      
      # Add a color bar to the legend
      legend_color_bar <- seq(min(heatmap_data), max(heatmap_data), length.out = 100)
      image.plot(x = rep(1, 100), y = legend_color_bar, z = matrix(legend_color_bar, nrow = 1), 
                 col = heat.colors(100), axes = FALSE, legend.only = TRUE, legend.width = 1.5)
      
      # Combine the heatmap and color bar
      layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1.5), heights = c(1, 4))
      return(heatmap_obj)
    })
    
    
    output$PCA_tab2 <- renderPlot({
      req(raw_data())
      
      timepoints <- sapply(colnames(raw_data()[-1]), timepoint_from_sample)
      
      replicates <- sapply(colnames(raw_data()[-1]), sample_replicate)
      #print(replicates)
      meta <- tibble(sample = colnames(raw_data()[-1]), timepoint = timepoints, replicate = replicates)
      
      dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(raw_data()[-1]),
        colData = as.data.frame(meta),
        design = ~1
      )
      
      pca_res <- prcomp(t(assay(dds)))
      pca_data <- as.data.frame(pca_res$x)
      pca_data$timepoint <- meta$timepoint
      
      pc_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
      pc_var_labels <- paste0("PC", 1:ncol(pca_res$x), ": ", round(pc_var, 0), "% Variance")
      
      selected_pc_x <- paste0("PC", input$pc_selection1)
      selected_pc_y <- paste0("PC", input$pc_selection2)
      
      plot <- ggplot(pca_data, aes_string(x = selected_pc_x, y = selected_pc_y, color = "timepoint")) +
        geom_point() +
        labs(title = paste("PCA Plot: Principal Components", selected_pc_x, "vs", selected_pc_y),
             x = paste0(selected_pc_x, ": ", round(pc_var[as.numeric(input$pc_selection1)], 0), "% Variance"),
             y = paste0(selected_pc_y, ": ", round(pc_var[as.numeric(input$pc_selection2)], 0), "% Variance"))
      
      return(plot)
    })
    
  }
  
  #run DESEQ on counts data-- tab 3               
  observeEvent(input$run_DE_dds, {
    # Run DESeq when the button is clicked
    req(input$file1, input$file2)
    
    counts_data <- raw_data()
    meta_data <- read.csv(input$file2$datapath)    
    
    rowData <- counts_data["gene"]
    counts_data <- as.data.frame(dplyr::select(counts_data, -gene))
    rownames(counts_data) <- rowData$gene
    
    colData <- DataFrame(
      samplename = colnames(counts_data),
      timepoint = sapply(colnames(counts_data), timepoint_from_sample)
    )
    
    colData <- colData[colData$timepoint %in% c('vP0', 'vAd', 'vP4', 'vP7'), ]
    
    selected_samples <- colData$samplename
    counts_data <- counts_data[, selected_samples]
    
    se <- SummarizedExperiment(
      assays = list(counts = as.matrix(counts_data)),
      colData = colData,
    )
    
    metadata(se) <- list(model = "model")
    
    unique_genes <- rownames(se)[!duplicated(rownames(se))]
    se_unique <- se[unique_genes, ]
    
    dds <- DESeqDataSet(se_unique, design = ~ timepoint)
    dds <- DESeq(dds)
   
    results <- results(dds, contrast = c("timepoint", "vAd", "vP0"))
    results_df <- as.data.frame(results)
    
    results_df$Gene <- rownames(results)
    results_df <- dplyr::select(results_df, Gene, everything())    #remove rownames
    
    result_list$results_df <- results_df
    result_list$dds <- dds
    return(result_list)
  })
  
  #OUPUT DESEQ TABLE
  output$DESEQ_res_table <- renderDT({
    req(input$file1, input$file2, input$run_DE_dds)
    results_df <- result_list$results_df
    datatable(results_df, options = list(columnDefs = list(list(targets = 0, visible = FALSE))))
  })
  
  #HELPER FUNCTION FOR VOLCANO PLOT-- creates volc_data
  label_res <- function(x) {
    
    labeled <- as_tibble(rownames = "genes", x) %>%
      mutate(volc_plot_status = case_when(
        x$padj < input$padj_threshold & x$log2FoldChange > 0 ~ "UP",
        x$padj < input$padj_threshold & x$log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )) 
    
    labeled_with_genes <- labeled %>% mutate(genes = rownames(x))
    
    volc_data <- labeled_with_genes %>% 
      dplyr::select(genes, volc_plot_status, everything())
    
    return(volc_data)
  }
  
  #REACTIVE VAL FOR PADJ SLIDER
  observe({
    req(input$file1, input$file2, input$run_DE_dds, input$padj_threshold)
    labeled_results <- label_res(result_list$results_df)
    
    labeled_results_filtered <- labeled_results[labeled_results$padj < input$padj_threshold, ]
    labeled_results_filtered <- labeled_results_filtered[complete.cases(labeled_results_filtered), ]
    
    filtered_results(tibble(labeled_results_filtered))
  })
  
  #CREATE VOLCANO PLOT
  output$Volcano_plot <- renderPlot({
    req(input$file1, input$file2, input$run_DE_dds, input$padj_threshold)
    
    labeled_results <- isolate(filtered_results())
    
    theme = labeled_results$volc_plot_status
    
    volc_plot <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = theme)) +
      geom_point(size = 1, alpha = 0.5) +
      labs(x = "log2FoldChange",
           y = "-log10(padj)",
           title = "Volcano plot of DESeq2 differential expression results vP0 vs. vAd)"
      ) + 
      scale_color_manual(values = c("UP" = "cornflowerblue", "DOWN" = "coral", "NS" = "forestgreen")) +
      theme_light() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black")
    
    return(volc_plot)
  })
 
  #VOLCANO PLOT FROM R SHINY ASSIGNMENT
  volcano_plot_rshy <-function(dataf, x_name, y_name, slider, color1, color2) {
    actual_data <- result_list$results_df
    actual_data <- transform(actual_data, padj = as.numeric(padj))
    
    ggplot(actual_data, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = !!sym(y_name) <= 10^slider)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(name = sprintf("-log10(%s) < %g", y_name, 10^slider), values = c(`TRUE` = color2, `FALSE` = color1)) +
      labs(
        title = sprintf("%s vs -log10(%s)", x_name, y_name),
        x = x_name,
        y = sprintf("-log10(%s)", y_name)
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",  
        legend.box = "horizontal",
        plot.title = element_text(face = "bold")
      )
  }

  #below observe blocks for R shiny plot
  observeEvent(input$go, {
    go_button_clicked(TRUE)
    data <- result_list$results_df
    last_plot(volcano_plot_rshy(data, input$button_x, input$button_y, input$padj_magnitude, input$color1, input$color2))
  })
  
  observeEvent({input$button_x; input$button_y; input$padj_magnitude}, {
    go_button_clicked(FALSE)
  })
  
  observe({
    req(result_list$results_df, go_button_clicked())
    if (go_button_clicked()) {
      data <- result_list$results_df
      last_plot(volcano_plot_rshy(data, input$button_x, input$button_y, input$padj_magnitude, input$color1, input$color2))
    }
  })
  
  output$rshiny_plot <- renderPlot({
    req(last_plot())
    last_plot()
  })
  
  observe({
    updateRadioButtons(session, "button_x", selected = "log2FoldChange")
    updateRadioButtons(session, "button_y", selected = "padj")
  })
  
  #create metadata for last tab
  meta_info_from_labels <- function(sample_names) {
    timepoints <- sapply(sample_names, timepoint_from_sample)
    replicates <- sapply(sample_names, sample_replicate)
    tble <- data.frame(sample = sample_names, timepoint = timepoints, replicate = replicates)    #this was a tibble
    return(tble)
  }
  
  
  #output the choice plot for last tab
  output$choice_plot <- renderPlot({
    req(input$go_choice, input$gene_search, input$meta_choice, input$file1)
    
    # Reactive val for gene names
    unique_genes <- rownames(norm_counts())
    
    # Reactive val for selected gene
    selected_gene <- grep(input$gene_search, unique_genes, value = TRUE, ignore.case = TRUE)
    
    samplenames <- colnames(norm_counts()[-1])
    meta_info_table <- meta_info_from_labels(samplenames)
    meta_info_table[, c("sample", "timepoint", "replicate")] <- lapply(meta_info_table[, c("sample", "timepoint", "replicate")], factor)
    
    trans_norm_counts <- t(norm_counts()[-1])
    
    # Merge metadata with transposed normalized counts data
    merged_data <- cbind(meta_info_table, trans_norm_counts)
    
    meta_choice_index <- which(colnames(merged_data) == input$meta_choice)
    
    # Filter data for the selected gene
    gene_data <- merged_data[, c(meta_choice_index, grep(selected_gene, colnames(merged_data)))]
    
    if (input$plot_choice == "box plot") {
      # Box plot
      p <- ggplot(gene_data, aes(x = as.factor(gene_data[, input$meta_choice]), y = gene_data[, selected_gene], fill = as.factor(gene_data[, input$meta_choice]))) +
        geom_boxplot() +
        labs(title = input$gene_search, x = input$meta_choice, y = "Normalized Counts") +
        theme(legend.position = "none")  # Optional: Hide legend if not needed
    } else if (input$plot_choice == "violin plot") {
      # Violin plot
      p <- ggplot(gene_data, aes(x = as.factor(gene_data[, input$meta_choice]), y = gene_data[, selected_gene], fill = as.factor(gene_data[, input$meta_choice]))) +
        geom_violin() +
        labs(title = input$gene_search, x = input$meta_choice, y = "Normalized Counts") +
        theme(legend.position = "none")  # Optional: Hide legend if not needed
    } else if (input$plot_choice == "beeswarm plot") {
      # Beeswarm plot
      p <- ggplot(gene_data, aes(x = as.factor(gene_data[, input$meta_choice]), y = gene_data[, selected_gene], color = as.factor(gene_data[, input$meta_choice]))) +
        geom_beeswarm() +
        labs(title = input$gene_search, x = input$meta_choice, y = "Normalized Counts") +
        theme(legend.position = "none")  # Optional: Hide legend if not needed
    } else if (input$plot_choice == "bar plot") {
      # Bar plot
      p <- ggplot(gene_data, aes(x = as.factor(gene_data[, input$meta_choice]), y = gene_data[, selected_gene], fill = as.factor(gene_data[, input$meta_choice]))) +
        geom_bar(stat = "identity") +
        labs(title = input$gene_search, x = input$meta_choice, y = "Normalized Counts") +
        theme(legend.position = "none")  # Optional: Hide legend if not needed
    }
    
    return(p)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
