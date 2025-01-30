library(shiny)
library(ggplot2)
library(shinyjs)
library(tidyr)
library(DT)
library(magrittr)
library(cowplot)
library(circlize)
library(RColorBrewer)
library(igraph)
library(plotly)
library(htmlwidgets)
library(htmltools)
library(tidyverse)
library(crosstalk)
library(dplyr)
library(scales)

#we are initializing the server from crosstalkserver initialize.r file


# Define server
server <- function(input, output, session) {
  
 

#####################################HOMEPAGE###########################################################
  # circos plot for homepage using parameter tuning cfunction
  output$circos_plot <- renderPlot({
    homepage_data <- parametertuning(realdata, input) # Call the parametertuning function here
    par(cex = 2)
    chordDiagram(homepage_data, transparency = 0.25, directional = FALSE, grid.col = category_color, annotationTrackHeight = c(0.05, 0.05),link.lwd = 1,    # Line width
                 link.lty = 1,    # Line type
                 link.border = 1)
  })
  
  # circo plot for emitter&reciever in homepage using calculate emittermatrix function
  output$circos_emitter <- renderPlot({
    data_filtered <- control_homepage(realdata, input)  # Call the control_homepage function
    homepage_emitter <- data_filtered %>%
      calculateEmitterMatrix()
    filteredrow <- SharedData$new(homepage_emitter)

    # Filter rows and columns with "U_" based on the checkbox 
    #check which divides the cell with U_ and regular 
    if (!is.null(input$filter_checkbox) && input$filter_checkbox) {
      filteredrow <- homepage_emitter[grepl("^U_", rownames(homepage_emitter)) & grepl("^U_", colnames(homepage_emitter)), grepl("^U_", colnames(homepage_emitter))]
    } else {
      filteredrow <- homepage_emitter[!grepl("^U_", rownames(homepage_emitter)) & !grepl("^U_", colnames(homepage_emitter)), !grepl("^U_", colnames(homepage_emitter))]
    }

    par(cex = 2)
    # Plot the chord diagram
    chordDiagram(
      filteredrow,
      transparency = 0.25,
      directional = TRUE,
      annotationTrackHeight = c(0.05, 0.05),
      link.lwd = 1,
      link.lty = 1,
      link.border = 1,
      direction.type = c("diffHeight", "arrows"),
      grid.col = category_color
    )
    
    circos.track(track.index = 2, panel.fun = function(x, y) {
      sector.index = get.cell.meta.data("sector.index")
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, niceFacing = TRUE)
    }, bg.border = NA)

    highlight.sector(rownames(filteredrow), track.index = 1, col = "red", 
                     text = "EMITTER", cex = 1.2, text.col = "white", niceFacing = TRUE,
                     values = 2,
                     sector.color = "red", sector.linewidth = 2,
                     xlim = c(0.1, 1), 
                     ylim = c(0, -100))
  })


  
  ############################################GENE INFORMATION#########################################
  #here we pull the gene name from gene_dictionary but filter the table using realdata
  # Filter the dataset based on gene search
  filtered_data <- reactive({
    req(input$gene_search)
 
    subset(gene_dictionary, grepl(input$gene_search, gene_name, ignore.case = TRUE)) %>%
      select(id_cp_interaction)
  })
  
  observeEvent(input$gene_search, {
    if (!is.null(input$gene_search)) {
      updateTabsetPanel(session, "tabs", selected = "Gene Information")
   }
  })
  
  
  

  
  # Render filtered data in the gene information table
  output$gene_table <- renderDT({
    req(filtered_data())
    
    # Filter inter_dic based on selected interactions #pulling the information for genetable from interdic
    selected_interactions <- filtered_data()$id_cp_interaction
    filtered_interactions <- inter_dic[inter_dic$id_cp_interaction %in% selected_interactions, c("id_cp_interaction", "ligand", "receptor", "ligand_pair")]
    
    datatable(
      filtered_interactions,
      options = list(
        pageLength = 10,
        rowCallback = JS(
          "function(row, data) {",
          "  $(row).on('click', function() {",
          "    var interactionID = data[0];",
          "    Shiny.setInputValue('interaction_select', interactionID);",
          "    Shiny.setInputValue('gene_table_row_clicked', new Date().getTime());",
          "    $('a[data-value=\"Interactions\"]').tab('show');",
          "  });",
          "}"
        )
      ),
      rownames = FALSE,
      escape = FALSE,
      selection = 'none'
      
    )
  })
  
  # Calculate cell_counts for all interactions in agene
  cell_counts_all <- reactive({
    req(filtered_data())
    selected_gene <- input$gene_search
    calculateCellCounts(selected_gene, input$cells_select, input$handlee, input$datasetss)
  })
  
   
 
  
  # Plot bar plot for all interactions in a gene
  output$analytics_plot_all <- renderPlot({
    cell_counts_for_plot <- cell_counts_all() %>%
      arrange(desc(n))
    
    cell_counts_for_plot$cell_ordered <- factor(cell_counts_for_plot$cell, levels = cell_counts_for_plot$cell)
    
    ggplot(cell_counts_for_plot, aes(x = cell_ordered, y = n, fill = cell_ordered)) +
      geom_bar(stat = "identity", color = "black", position = "identity", width = 0.8) +
      scale_fill_manual(values = category_color) +
      labs(x = "Cell Type", y = expression(bold("Frequency")), title = NULL, fill = "Cell") +
      ggtitle("Interactions for the Gene") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
      geom_text(aes(label = n, fontface = "bold"), vjust = -0.5) + scale_x_discrete(drop = FALSE)
  })
  


  # Calculate adjacency table for circoplot in genetab
  gene_pair_counts <- reactive({
    selected_gene <- input$gene_search
    
    data <- realdata %>%
      filter(gene_a == selected_gene | gene_b == selected_gene) %>%
      control_genepage(input)
    
    gene_counts <- calculateadjacencymatrix(data)
    
    # Calculate row sums
    row_sums <- rowSums(gene_counts)
    
    # Sort rows based on row sums
    gene_counts <- gene_counts[order(row_sums, decreasing = TRUE), ]
    #calculating for circoplot on genepage using function calculateemittermatrix
    emitter_genepage <- realdata %>%
      filter(gene_a == selected_gene | gene_b == selected_gene) %>%
      control_genepage(input) 
    
    emitter_genepage <- calculateEmitterMatrix(emitter_genepage)
    
    
    # Get the row and column names of the matrix
    row_names <- rownames(emitter_genepage)
    col_names <- colnames(emitter_genepage)
    
    # Use grepl on row and column names
    filtered_rows <- grepl("^U_", row_names)
    filtered_cols <- grepl("^U_", col_names)
    
    # Subset the matrix using the logical vectors
    filtered_emitter_genepage <- emitter_genepage[!filtered_rows, !filtered_cols]
    
    list(gene_counts = gene_counts, emitter_genepage = filtered_emitter_genepage)
  })
  
  
  # plot_chordDiagram <- function(count_mat, category_color, symmetric = T) {
  #   chordDiagram(count_mat, transparency = 0.4, grid.col = category_color, link.arr.width = 8, annotationTrackHeight = c(0.07, 0.07))
  # }
  # 
  
  # Render circos plot in the Analytics tab
  output$circos_plot_interactions <- renderPlot({
    par(cex = 1.4)
    
    if (input$Gene_emitter_checkbox) {
      chordDiagram(gene_pair_counts()$emitter_genepage, transparency = 0.4, grid.col = category_color, link.arr.width = 8, annotationTrackHeight = c(0.07, 0.07))
      highlight.sector(rownames(gene_pair_counts()$emitter_genepage), track.index = 1, col = "red", 
                       text = "EMITTER", cex = 1.2, text.col = "white", niceFacing = TRUE,
                       values = 2,
                       sector.color = "red", sector.linewidth = 2,
                       xlim = c(0.1, 1), 
                       ylim = c(0, -100))
      
    } else {
      chordDiagram(gene_pair_counts()$gene_counts, transparency = 0.4, grid.col = category_color, link.arr.width = 8, annotationTrackHeight = c(0.07, 0.07))
    }
  })
  
 

  
  
  # Render DataTable for gene_pair_counts
  output$gene_pair_counts <- renderDT({
    if (input$Gene_emitter_checkbox) {
      datatable(gene_pair_counts()$emitter_genepage, rownames = TRUE, options = list(dom = 't'), selection = 'none')
    } else {
      datatable(gene_pair_counts()$gene_counts, rownames = TRUE, options = list(dom = 't'), selection = 'none')
    }
  })
  
  
  

  # Render the selected interaction details
  output$selected_interaction <- renderText({
    selected_interaction <- input$interaction_select
    searched_gene <- input$gene_search
    if (!is.null(selected_interaction)) {
      matching_ligand_pair <- unique(inter_dic$ligand_pair[inter_dic$id_cp_interaction == selected_interaction])
      paste0(selected_interaction, " : ",matching_ligand_pair)
    }
  })
  
  # Function to create hyperlink URL
  createPubMedLink <- function(url) {
    pubmed_id <- gsub("[^0-9]", "", url)
    if (nchar(pubmed_id) > 0) {
      hyperlink <- sprintf('<a href="https://pubmed.ncbi.nlm.nih.gov/%s" target="_blank">%s</a>', pubmed_id, url)
    } else {
      hyperlink <- url
    }
    hyperlink
  }
  
  output$interaction_details_table <- renderDataTable({
    selected_interaction <- input$interaction_select
    
    if (!is.null(selected_interaction)) {
      interaction_data <- inter_dic[inter_dic$id_cp_interaction == selected_interaction, ]
      
      data <- data.frame(
        
        Ligand = unique(inter_dic$ligand[inter_dic$id_cp_interaction == selected_interaction]),
        Receptor = unique(inter_dic$receptor[inter_dic$id_cp_interaction == selected_interaction]),
        Curated = ifelse(unique(inter_dic$annotation_strategy[inter_dic$id_cp_interaction == selected_interaction]) == "curated", "Yes", "No"),
        Reference = unique(inter_dic$source[inter_dic$id_cp_interaction == selected_interaction])
      )
      
      # Apply formatting to the Reference column to create clickable hyperlinks
      data$Reference <- sapply(data$Reference, createPubMedLink)
      
      
      
      
      datatable(data, rownames = FALSE, escape = FALSE,
                options = list(dom = 't', columnDefs = list(list(targets = "_all", className = "dt-center"))),selection = 'none',
                callback = JS("table.column(3).nodes().to$().css({cursor: 'pointer', color: 'blue', 'text-decoration': 'underline'});"),)
    }
  })
  
  


  
  
                                                                 ######INTERACTIONS TAB#######
  
    
  # Calculate cell counts and plot graph for the selected interaction ID
  observe({
    input$gene_table_row_clicked
    
    req(input$interaction_select)
    # Define the list of variables
    desired_cells <- c("Astrocytes", "Microglia", "NeuronE", "NeuronI", "OPC", "Oligos", "Endothelial")
    # Calculate cell counts
    cell_counts <- realdata %>%
      filter(id_cp_interaction == input$interaction_select) %>%
      control_interactionspage(input) # Assign the result to a variablem
    
    cell_counts <- cell_counts %>%
      select(id_cp_interaction, cell1, cell2) %>%
      mutate(row_index = 1:n()) %>% 
      gather("foo", "cell", cell1, cell2) %>%
      select(-foo) %>% 
      distinct() %>% 
      count(cell, id_cp_interaction) %>%
      complete(cell = desired_cells, id_cp_interaction, fill = list(n = 0)) %>%
      left_join(data.frame(cell = desired_cells), by = "cell") %>%
      arrange(id_cp_interaction, cell)
    
    analytics_emitter <- realdata %>%
      filter(id_cp_interaction == input$interaction_select) %>%
      control_interactionspage(input) # Assign the result to a variable
    
   
    #emitter matrix
    analytics_emitter <- calculateEmitterMatrix(analytics_emitter)
    
    # Plot graph
    cell_counts_for_plot <- cell_counts %>%
      arrange(desc(n))
    
    cell_counts_for_plot$cell_ordered <- factor(cell_counts_for_plot$cell, levels = cell_counts_for_plot$cell)
    
    output$analytics_plot <- renderPlot({
      selected_gene <- unique(filtered_data()$interacting_pair)
      gene_name <- unique(filtered_data()$gene_name)
      
     
      
      
      ggplot(cell_counts_for_plot, aes(x = cell_ordered, y = n, fill = cell_ordered)) +
        geom_bar(stat = "identity", color = "black", position = "identity", width = 0.8) +
        scale_fill_manual(values = category_color) +
        labs(x = "Cell Type", y = expression(bold("Frequency")), title = NULL, fill = "Cell") +
        ggtitle(bquote(atop(bold(.(input$interaction_select)), atop(bold(.(gene_name)), italic(.(selected_gene)))))) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
              axis.text.y = element_text(size = 12,face="bold"),
              axis.title.y = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14,face = "bold"),
              legend.title = element_text(size = 16,face = "bold"),
              plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
        geom_text(aes(label = n, fontface = "bold"), vjust = -0.5)
    })
    
    
    
    # Calculate adjacency table
    cell_pair_counts <- realdata %>%
      filter(id_cp_interaction == input$interaction_select) %>% 
    control_interactionspage(input)
    
    cell_pair_counts <- calculateadjacencymatrix(cell_pair_counts)
    
    # Calculate row sums
    row_sums <- rowSums(cell_pair_counts)
    
    # Sort rows based on row sums
    cell_pair_counts <- cell_pair_counts[order(row_sums, decreasing = TRUE), ]
    
    output$circos_plot_analytics <- renderPlot({
      par(cex = 1.4)
      if (input$analytics_emitter_checkbox) {
        chordDiagram(analytics_emitter, transparency = 0.4, grid.col = category_color, link.arr.width = 8, annotationTrackHeight = c(0.07, 0.07))
      } else {
        chordDiagram(cell_pair_counts, transparency = 0.4, grid.col = category_color, link.arr.width = 8, annotationTrackHeight = c(0.07, 0.07))
      }
    })
    
    # Create the density plot
    output$density_plot <- renderPlot({


      # Choose a specific handle and interaction ID
      selected_interaction <- input$interaction_select
      selected_handle <- input$handles

      # Filter the gene_dictionary dataframe to retrieve gene names associated with the selected interaction
      selected_gene_names <- gene_dictionary %>%
        filter(id_cp_interaction == selected_interaction) %>%
        pull(gene_name)

      # Now you can use the selected_gene_names in your code
      # Filter the gene expression data for the selected handle and the extracted gene names
      selected_handle_expr <- gex_list[[selected_handle]]
      selected_gene_expr <- selected_handle_expr %>%
        filter(gene %in% selected_gene_names)

       if (input$expr_percentiles) {
        # Create a percentile plot
         ggplot(selected_handle_expr, aes(x = percentile)) +
           geom_rect(aes(xmin = -Inf, xmax = Inf, ymin =0.3, ymax = 0.7), alpha = 0.2) +
           geom_vline(aes(xintercept = percentile, color = gene), data = selected_gene_expr, linewidth=1.5) +
           labs(
             x = "Percentile",
             y = NULL,  # No y-axis label
             title = "Gene Expression Percentiles by Cell Type"
           ) + scale_y_continuous(limits = c(0,1) , breaks = NULL) +
           facet_wrap(~ celltype, scales = "free")+
           theme_minimal() + theme(panel.grid = element_blank()) +
           theme(
             legend.text = element_text(size = 12, face = "bold"),
             legend.title = element_text(size = 14, face = "bold") ,
             axis.title.x = element_text(size = 20, face = "bold")
           ) +
           theme(strip.text = element_text(size=9 , face = "bold"))
      } else {
        # Create a density plot of gene expression
        ggplot(selected_handle_expr, aes(x = counts ,fill=celltype)) +
          geom_density(alpha = 0.6) +
          geom_vline(aes(xintercept = counts, color = gene), data = selected_gene_expr, linewidth=1.5) +
          labs(x = "Gene Expression", y = "Density" ,title = "Gene Expression Density Plot for Brase et.al ") +facet_grid(~ celltype, scales = "free", switch = "x")+
          facet_wrap(~ celltype, ncol = 1, scales = "free") +
          scale_x_continuous(trans = 'log10' ,labels = comma) +
          theme_minimal() +
          theme(panel.grid = element_blank()) +
          guides(fill = guide_none()) +
          theme(
            legend.text = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 14, face = "bold") ,
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold")
          ) +
          theme(strip.text = element_text(size=9 , face = "bold"))
      }
    })

   
    
   
    
    # Render DataTable for cell_pair_counts
    output$cell_pair_counts <- renderDT({
      datatable(cell_pair_counts, rownames = TRUE, options = list(dom = 't'),selection ='none' )
    })
    
    
  })
  
  # Reload the app when coming back to the "Gene Information" tab
  observeEvent(input$gene_information_tab, {
    session$reload()
  })
  
}