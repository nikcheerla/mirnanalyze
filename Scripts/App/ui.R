library(shiny)
library(stats)

# Define the overall UI
shinyUI(
  
  # Use a fluid Bootstrap layout
  fluidPage(    
    
    # Give the page a title
    titlePanel("microRNA based Cancer Diagnosis and Treatment Recommendation"),
    
    # Generate a row with a sidebar
    sidebarLayout(      
      
      # Define the sidebar with one input
      sidebarPanel(
        
        img(src = "mirna.jpg", height = 200, width = 200),
        
        h4("Inputs"),
        p("In this application, enter miRNA data in the csv file format along with some clinical information to get diagnosis"),
        textInput("Age", label = h5("Age of the Patient"), 
                  value = "50"),
        
        
        selectInput("Gender", 
                    label = h5("Gender of the Patient"),
                    choices = c("Male","Female"),
                    selected = "Male"),
        
        selectInput("Ethnicity", 
                    label = h5("Ethnicity of the Patient"),
                    choices = c("White","Latino", "African American", "American Indian ", "Asian", "Other"),
                    selected = "Other"),
        
        p("The first column of this CSV file should have miRNA name in hsa- format. Second column contains miRNA values"),
        
        
        tags$hr(),
        checkboxInput('header', 'Header', TRUE),
        checkboxInput('log2', 'Log2 Normalized?', TRUE),
        radioButtons('sep', 'Separator',
                     c(Comma=',',
                       Semicolon=';',
                       Space = ' ',
                       Tab='\t'),
                     ','),
        
        fileInput('file1', 'Choose CSV File',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.CSV',
                           '.txt',
                           '.TXT',
                           '.csv')), 
        
        checkboxInput('upload', 'Check to give us permission to upload this data to our database', FALSE)
        
        
      ),
      
      
      # Create a spot for the barplot
      mainPanel(
        h4("About This Tool"), 
        p("This web tool lets the user upload patient 
          information including clinical data and micro RNA data. Then, 
          using built-in predictive models, 
          it diagnoses the patient as well as recommends the best possible 
          treatment with remission probability."),
        
        p("The tool can currently diagnose between 21 different kinds of cancers and recommend treatments for 9 different kinds of cancers. If the diagnosis is for a cancer for which there is no treatment recommendation, the tool will notify that case."),
        
        p("After you upload the miRNA expression and the clinical data, wait 1-2 minutes while the computation occurs."),
        br(),
        
        h4("Diagnosis and Treatment Recommendation"),
        p(" This may take a few minutes. Please wait..."),
        
        textOutput("text1"),
        br(),
        htmlOutput("text2"),
        
        br(),
        h4("A few entries of the microRNA data you uploaded:"),
        tableOutput('contents1'), 
        
        br(),
        h4("Log2 Normalized Data Before Imputation"),
        p("Note that the tool only extracts a subset of 60 microRNAs from the used uploaded data. It less than 20% of those miRNA has unknown values, the tool imputes those values and can proceed further. If more than 20% of miRNA are missing, the tool will stop further processing"),
        
        tableOutput('contents'), 
        br()
        
     
        
        )
      
      )
    )
  )