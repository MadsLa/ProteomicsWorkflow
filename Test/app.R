# Load necessary libraries
library(shiny)
library(quarto)

# Define UI
ui <- fluidPage(
  titlePanel("R Shiny Quarto PDF Report Generation Example"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("input_file1", "Choose File 1", accept = c(".csv", ".txt")),
      fileInput("input_file2", "Choose File 2", accept = c(".csv", ".txt")),
      textInput("text_1", "Enter Text 1", ""),
      textInput("text_2", "Enter Text 2", ""),
      actionButton("generate_report", "Generate PDF Report")
    ),
    
    mainPanel(
      h4("Uploaded Files"),
      textOutput("file1_name"),
      textOutput("file2_name"),
      
      h4("Entered Text"),
      textOutput("text1_output"),
      textOutput("text2_output"),
      
      downloadButton("download_report", "Download PDF Report")
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$file1_name <- renderText({
    if (is.null(input$input_file1)) {
      "No file uploaded for File 1"
    } else {
      paste("File 1 uploaded:", input$input_file1$name)
    }
  })
  
  output$file2_name <- renderText({
    if (is.null(input$input_file2)) {
      "No file uploaded for File 2"
    } else {
      paste("File 2 uploaded:", input$input_file2$name)
    }
  })
  
  output$text1_output <- renderText({
    paste("Text 1:", input$text_1)
  })
  
  output$text2_output <- renderText({
    paste("Text 2:", input$text_2)
  })
  
  # Create a reactive value to store the report path
  report_path <- reactiveVal()
  
  observeEvent(input$generate_report, {
    # Create a temporary directory to store files
    temp_dir <- tempdir()
    
    # Copy the Quarto template file to the temporary directory
    file.copy("template.qmd", file.path(temp_dir, "template.qmd"))
    
    # Define file paths
    input_file1_name <- if (is.null(input$input_file1)) "None" else input$input_file1$name
    input_file2_name <- if (is.null(input$input_file2)) "None" else input$input_file2$name
    
    # Create a temporary R script to render the Quarto document
    r_script_content <- paste0(
      "input_file1 <- '", input_file1_name, "'\n",
      "input_file2 <- '", input_file2_name, "'\n",
      "text_1 <- '", input$text_1, "'\n",
      "text_2 <- '", input$text_2, "'\n",
      "quarto::quarto_render('", file.path(temp_dir, "template.qmd"), "', output_format = 'pdf')"
    )
    
    # Save the R script to a temporary file
    r_script_path <- file.path(temp_dir, "render_report.R")
    writeLines(r_script_content, con = r_script_path)
    
    # Run the R script to create the PDF
    source(r_script_path)
    
    # Define PDF file path
    pdf_file <- file.path(temp_dir, "template.pdf")
    
    # Store the path of the PDF file
    report_path(pdf_file)
  })
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste("report-", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      file.rename(report_path(), file)
    },
    contentType = "application/pdf"
  )
}

# Run the application 
shinyApp(ui = ui, server = server)