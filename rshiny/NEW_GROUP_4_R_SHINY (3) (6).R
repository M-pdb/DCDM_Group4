#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Install required packages if they are not already installed
if (!requireNamespace("RMySQL", quietly = TRUE)) install.packages("RMySQL")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("dplyr)", quietly = TRUE)) install.packages("shiny")

#Load the required packages
library(shiny)
library(RMySQL)
library(dplyr)
library(DBI)
library(ggplot2)
library(tidyr)
library(dplyr)

#Defines the useer interface
ui <- fluidPage( #creates the main container for the UI
  titlePanel("Mouse Phenotype Analysis"), #Sets the title of the application
  tabsetPanel( #Creates the actual tabbed interface
    tabPanel("Knockout Mouse Analysis", #Sets the title of the first plot
             sidebarLayout( #Divides the tab into a sidebar and the main area
               sidebarPanel( #Creates the dropdown menu
                 selectInput("knockoutMouse", "Select Knockout Mouse:", choices = c("Loading..." = ""))
               ),
               mainPanel( #This creates the main plot
                 plotOutput("knockoutPlot", height = "800px") #this sets the main plot
               )
             )
    ),
    tabPanel("Phenotype Group Analysis", #This creates the second tab, for the second plot, called Phenotype Group Analysis
             sidebarLayout(
               sidebarPanel(
                 selectInput("phenoGroup", "Select Phenotype Group:", choices = NULL) #This creates the drop down menu, containing the groups
               ),
               mainPanel(
                 plotOutput("phenoGroupPlot")
               )
             )
    ),
    tabPanel("Heatmap", 
             plotOutput("heatmapPlot", height = "600px", width = "200%") #This creates the panel for the heatmap, this is not interactive, the width is set at 200% so that more gene symbols can be added
    )
  )
)

server <- function(input, output, session) {
  # Database connection function with error handling
  get_data <- function(query) {
    tryCatch({
      con <- dbConnect(RMySQL::MySQL(),
                       dbname = "mouse_genomev3", #This depends on the name of the database
                       host = "localhost",
                       port = 3306,
                       user = "root", #This is the user name, which is usually set to default, which is root
                       password = "password") # this is where the password for SQL needs to be put in
      on.exit(dbDisconnect(con))
      results <- dbGetQuery(con, query)
      return(results)
    }, error = function(e) { 
      message("Database error: ", e$message) #if there is a problem with the connection, this will say that there is an error
      return(NULL) 
    })
  }
  
  # Get knockout mouse options
  knockout_options <- reactive({
    query <- " #this preforms a SQL query and stores it in query
      SELECT DISTINCT ad.mouse_strain, ad.gene_symbol, g.gene_accession_id #This uses a SQL query, which selects unique mouse strains, gene symbols and gene acession_ids 
      FROM Analysis_Data ad 
      JOIN Genes g ON ad.gene_accession_id = g.gene_accession_id #makes an inner join between the Genes (g) table and the Aanalysis_Data (ad). gene accession_id is the link between the two tables
      ORDER BY ad.mouse_strain, ad.gene_symbol #the results are sorted by mouse strain in the Analysis_Data and then be gene_symbol, which is also found in the Analysis_Data
    "
    results <- get_data(query) #this gets the data from the query, and stores it in results
    if (is.null(results) || nrow(results) == 0) { #Error handling step. Checks if there is NULL which would mean that the query has failed, or if there are zero rows, meaning that there is no data
      return(character(0))
    }
    return(setNames( #This creates a vector for the use in the uswer interface
      paste(results$mouse_strain, results$gene_symbol, results$gene_accession_id, sep = "|"),
      paste("Mouse Strain: (", results$mouse_strain, ") & Knocked Out Gene(", results$gene_symbol, ")", sep = "")
    ))
  })
  
  # Update knockout mouse choices
  observe({
    choices <- knockout_options()
    if (length(choices) > 0) {
      updateSelectInput(session, "knockoutMouse", choices = choices)
    } else {
      updateSelectInput(session, "knockoutMouse", choices = c("No options available" = ""))
    }
  })
  
  # Reactive data for knockout mouse analysis
  knockout_data <- reactive({
    req(input$knockoutMouse)
    selected <- strsplit(input$knockoutMouse, "\\|")[[1]]
    mouse_strain <- selected[1]
    gene_symbol <- selected[2]
    gene_accession_id <- selected[3]
    query <- sprintf("
      SELECT gp.group_name, ad.parameter_name, ad.pvalue
      FROM Group_Parameters gp
      JOIN IMPC_Parameter_Description ipd ON gp.parameter_id = ipd.parameter_id #Preforms a inner join between IMPC_Parameter_Description table (ipd) and the Group_Parameter table (gp). The parameter_id is the link between the two tables
      JOIN Analysis_Data ad ON ipd.parameter_id = ad.parameter_id #Makes an inner join between the Analysis_Data table (ad) and the IMPC_Parameter_Description table (ipd). The parameter_id is the link between the two tables
      WHERE ad.mouse_strain = '%s' AND ad.gene_symbol = '%s' AND ad.gene_accession_id = '%s' #Filters for mouse strains, gene_symbols and gene accesission ids
      ORDER BY gp.group_name, ad.pvalue #sorted alphabetically by group name, and then within each group the results are further sorted by their p values
    ", mouse_strain, gene_symbol, gene_accession_id)
    data <- get_data(query)
    if (!is.null(data) && nrow(data) > 0) {
      data <- data %>%
        filter(!grepl("uncategorized|uncharacterized", tolower(group_name)))
      data$neg_log10_p <- -log10(data$pvalue)
      return(data)
    } else {
      return(NULL)
    }
  })
  
  # Reactive data for phenotype group analysis
  pheno_group_data <- reactive({
    query <- "
      SELECT ad.gene_symbol, ad.gene_accession_id, ad.pvalue, gp.group_name, ad.parameter_name
      FROM Analysis_Data ad
      JOIN IMPC_Parameter_Description ipd ON ad.parameter_id = ipd.parameter_id
      JOIN Group_Parameters gp ON ipd.parameter_id = gp.parameter_id
      WHERE ad.gene_symbol COLLATE utf8mb4_bin IN ('Itln1', 'Lamb3', 'Mapkap1', 'Ttc39b')
        AND gp.group_name NOT LIKE '%Uncategorized%'
      ORDER BY ad.gene_symbol, gp.group_name, ad.pvalue;
    "
    data <- get_data(query)
    data$neg_log10_p <- -log10(data$pvalue) #Calculates the negative log base of the p value
    return(data)
  })
  
  # Update phenotype group choices
  observe({ #This function creates a reactive observer, this will re run whenever the input changes
    pheno_groups <- unique(pheno_group_data()$group_name)
    updateSelectInput(session, "phenoGroup", choices = pheno_groups, selected = pheno_groups[1])
  })
  
  # Knockout Mouse Analysis Plot. #This creates the first plot, which is a mixture of a box plot and a jitter plot
  output$knockoutPlot <- renderPlot({
    data <- knockout_data()
    req(data)
    ggplot(data, aes(x = group_name, y = neg_log10_p, color = group_name)) +
      geom_jitter(width = 0.3, alpha = 0.5) + #adds the jitter plot
      geom_boxplot(alpha = 0.3, outlier.shape = NA) + #adds the box plot
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") + #No legend was added to make it look cleaner
      labs(title = paste("Phenotypes for", gsub("\\|.*$", "", input$knockoutMouse)), #Adds labels to the plot
           x = "Phenotype Groups",
           y = "-log10(p-value)") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      annotate("text", x = 1, y = -log10(0.05), label = "p = 0.05", vjust = -0.5, color = "red")
  })
  
  # Phenotype Group Analysis Plot. This creates a box plot
  output$phenoGroupPlot <- renderPlot({
    req(input$phenoGroup)
    data <- pheno_group_data() %>%
      filter(group_name == input$phenoGroup)
    ggplot(data, aes(x = gene_symbol, y = neg_log10_p, fill = gene_symbol)) +
      geom_boxplot() + #creates the boxplot
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Analysis of", input$phenoGroup, "For Genotypes of Interest"), #Adds the labels 
           x = "Genotypes of Interest",
           y = "-log10(p-value)",
           fill = "Genotype of Interest") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      annotate("text", x = 1, y = -log10(0.05), label = "p = 0.05", vjust = -0.5, color = "red")
  })
  
  # This makes the Heatmap Plot for the third plot
  output$heatmapPlot <- renderPlot({
    query <- "
    SELECT ad.gene_symbol, ad.pvalue, gp.group_name 
    FROM Analysis_Data ad
    JOIN IMPC_Parameter_Description ipd ON ad.parameter_id = ipd.parameter_id
    JOIN Group_Parameters gp ON ipd.parameter_id = gp.parameter_id
    WHERE ad.gene_symbol COLLATE utf8mb4_bin IN ('Itln1', 'Lamb3', 'Mapkap1', 'Ttc39b');
    "
    data <- get_data(query)
    req(data)
    
    heatmap_data <- data %>%
      pivot_wider(names_from = gene_symbol, values_from = pvalue, values_fn = mean) #Turns the into long format, which is more suitable for a heatmap. And if there are multiple p values for a gene, the mean is
    
  
    par(mar = c(20, 4, 4, 20), xpd = TRUE)
    heatmap(as.matrix(heatmap_data[,-1]),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        labRow = heatmap_data$group_name,
        labCol = colnames(heatmap_data)[-1], # Ensure column labels are included
        main = "P-values of genes in knock_out Mice Heatmap",
        scale = "column",
        margins = c(10, 10))
    
    
  }, height = 700 , width = 700)
}

shinyApp(ui, server)
