library(DT)
library(shiny)
library(circlize)

ui <- fluidPage(
  mainPanel(
    tabsetPanel( 
      id = "tabs",
      tabPanel(
        "OVERVIEW",
        fluidRow(
          column(
            width = 6,
            sidebarPanel(
              selectInput("gene_search", "Gene", choices = unique(gene_dictionary$gene_name), multiple = TRUE)
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            div(
              style = "position: relative;",
              div(
                style = "position: absolute; top: 20px; right: -200px;",
                selectInput("cell_selection", "Select Cells", choices = c("all", "Astrocytes", "Microglia", "NeuronE", "NeuronI", "OPC", "Oligos", "Endothelial"))
              ),
              div(
                style = "position: absolute; top: 90px; right: -200px;",
                selectInput("handle", "Handle", choices = c("all", "CO", "AD", "ADAD", "TREM2"))
              ),
              div(
                style = "position: absolute; top: 160px; right: -200px;",
                selectInput("dataset", "Dataset", choices = c("Brase", "Colonna", "Lau", "Morabito", "all")),
                div(
                  style = "position: absolute; top: 90px; right: -220px; font-size: 14px; color: #888; border: 1px solid #ccc; padding: 10px; background-color: #f7f7f7; max-width: 600px;",
                  style = "font-size: 16px; font-weight: bold;",
                  "Handle",
                  tags$ul(
                    style = "list-style-type: none; padding-left: 0;",
                    tags$li(tags$strong(" CO - Control")),
                    tags$li(tags$strong(" AD - Alzheimer's disease")),
                    tags$li(tags$strong(" ADAD - Autosomal Dominant Alzheimer's Disease (ADAD)")),
                    tags$li(tags$strong(" TREM2 - Triggering Receptor Expressed On Myeloid Cells 2"))
                  )
                )
              ),
              plotOutput("circos_plot", height = "870px", width = "100%")
            )
          )
        ),
        fluidRow(
          column(
            width = 12,
            div(
              style = "position: relative;",
              checkboxInput("filter_checkbox", "Non Significant interactions"),
              plotOutput("circos_emitter", height = "800px", width = "100%")
            )
          )
        )
      ),
      tabPanel(
        "Gene Information",
        fluidRow(
          column(
            width = 12,
            DTOutput("gene_table"),
            column(
              width = 6,
              div(
                style = "position: relative;",
                plotOutput("analytics_plot_all"),
                div(
                  style = "position: absolute; top: 240px; right: -860px;",
                  selectInput("cells_select", "Select Cells", choices = c("all", "Astrocytes", "Microglia", "NeuronE", "NeuronI", "OPC", "Oligos", "Endothelial"))
                ),
                div(
                  style = "position: absolute; top: 90px; right: -860px;",
                  selectInput("handlee", "Handle", choices = c("CO", "all", "AD", "ADAD", "TREM2"))
                ),
                div(
                  style = "position: absolute; top: 160px; right: -860px;",
                  selectInput("datasetss", "Dataset", choices = c("Brase","all",  "Colonna", "Lau", "Morabito"))
                )
              )
            ),
            column(
              width = 6,
              div(
                style = "position: relative;",
                checkboxInput("Gene_emitter_checkbox", "Emitter & Reciever"),
                plotOutput("circos_plot_interactions", height = "600px", width = "100%")
              )
            )
          )
        ),
        
        DTOutput("gene_pair_counts")
        
      ),
      tabPanel(
        "Interactions",
        fluidRow(
          column(
            width = 12,
            h1(textOutput("selected_interaction"), style = "text-align: center; font-weight: bold;")
          ),
          column(
            width = 12,
            div(
              style = "position: relative;",
              dataTableOutput("interaction_details_table")
            )
          ),
          column(
            width = 6,
            div(
              style = "position: relative;",
              plotOutput("analytics_plot"),
              div(
                style = "position: absolute; top: 240px; right: -860px;",
                selectInput("cells_selections", "Select Cells", choices = c("all", "Astrocytes", "Microglia", "NeuronE", "NeuronI", "OPC", "Oligos", "Endothelial"))
              ),
              div(
                style = "position: absolute; top: 90px; right: -860px;",
                selectInput("handles", "Handle", choices = c("all", "CO", "AD", "ADAD", "TREM2"))
              ),
              div(
                style = "position: absolute; top: 160px; right: -860px;",
                selectInput("datasets", "Dataset", choices = c("all", "Brase", "Colonna", "Lau", "Morabito"))
              )
            )
          ),
          column(
            width = 6,
            div(
              style = "position: relative;",
              checkboxInput("analytics_emitter_checkbox", "Emitter & Reciever"),
              plotOutput("circos_plot_analytics", height = "600px", width = "100%")
            )
          )
        ),
        DTOutput("cell_pair_counts"),
        checkboxInput("expr_percentiles","Expression Percentiles"),
        plotOutput("density_plot", width = "100%", height = "1200px")
      ),
      tabPanel(
        "Description",
        HTML(
          '<p><b>Explore the intricate cellular interactions and gene expression patterns within the human brain with our dedicated platform. Delve into the world of cellular communication, neurodegenerative diseases like Alzheimer\'s, and the role of various brain cell types in maintaining brain health. Uncover insights through single-nucleus transcriptomic analysis, genetic research, and gene co-expression networks. Join us in unraveling the complexities of brain function and its connection to diseases through comprehensive research and exploration.</b></p>
     <p>Our website offers a comprehensive exploration of the following key aspects:</p>
     <ul>
     <br>
       <li><b>Cellular Crosstalk in the Brain:</b><br>
       We uncover the intricate network of communication established among different brain cell types. Through single-nucleus transcriptomic analysis, we reveal the ligand-receptor interactions that drive cellular communication, providing insights into how these interactions influence brain physiology and pathology.</li>
       <br>
       <li><b>Understanding Neurodegenerative Diseases:</b><br>
       We focus on Alzheimer\'s disease (AD) as a prime example of how aberrant cellular crosstalk can contribute to neurodegeneration. By analyzing genetic risk loci and functional genomics data, we unravel the role of microglia, the brain\'s resident immune cells, in mediating AD genetic risk. We explore how disrupted crosstalk between microglia and other cell types is linked to AD development.</li>
       <br>
       <li><b>Reconstructing Cellular Networks:</b><br>
       Our website presents a systematic approach to reconstructing the gene co-expression networks that underlie cellular crosstalk. We highlight how specific crosstalk interactions involving AD-related genes are enriched in microglia-neuron communication, shedding light on key players in AD pathology.</li>
       <br>
       <li><b>Microglia Activation and Disease Progression:</b><br>
       We investigate the role of critical crosstalk interactions, such as the SEMA6D-TREM2 axis, in modulating microglia activation. By analyzing gene expression patterns in relation to disease severity, we uncover how these interactions are disrupted in advanced AD stages.</li>
       <br>
       <li><b>Spatial Insights and In Vitro Validation:</b><br>
       Our platform leverages spatial transcriptomics data to provide insights into the expression patterns of key genes in proximity to Aβ plaques, a hallmark of AD. We validate our findings using in vitro human iPSC-derived microglia, demonstrating the functional impact of crosstalk interactions on microglia activity and cytokine release.</li>
     </ul>'
        )
      )
      
    )
  )
)