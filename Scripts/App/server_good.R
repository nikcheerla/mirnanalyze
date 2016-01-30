library(shiny)
library(stats)
library("e1071")
library("preprocessCore")
library("caTools")
library("Hmisc")
library("impute")
library("glmnet")
library("shiny")
library("rdrop2")
library("caret")
library("kernlab")
library("Matrix")
library("foreach")

load("mirnaData_reduced_60.RData")
load("clinicalData.RData")
load("clinFit60.RData")
load("svmFit60.Rdata")

source("helpers.R")

x <- mirna1
y <- class
token <- readRDS("droptoken.rds")



# Define a server for the Shiny app
shinyServer(function(input, output, session) {
  
  #################################################################################
  
  
  # read the miRNA csv file uploaded by the user
  dataInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    inData <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                       quote=input$quote)
    
    #return(inData)
    
  })
  
  # print out the first three rows of the input so users can check sanity
  output$contents1 <- renderTable({
    toplot <- dataInput()
    toplot <- toplot[1:3, ]
  })
  
  # output the pre processed data after imputation
  output$contents <- renderTable({
    
    indata <- dataCleanup();
    if (!is.null(indata)) {
      wantedData <- indata[[1]]
      print("here")
      colnames(wantedData) <- c("miRNA name", "miRNA Value")
      wantedData <- wantedData
    }
    
  })
  
  
  showProgress <- reactive({ 
    
    data <- dataInput();
    
    if (!is.null(data)) { 
      ################ progress bar ################################
      progress <- shiny::Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      
      progress$set(message = 'Pre Processing',
                   detail = 'This may take a while...')
      
      for (i in 1:50) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      ##############################################################
      
      
    }
    
  })
  
  # matches the input data to mirna1 matrix. 
  dataCleanup <- reactive({
    
    data <- dataInput()
    if (!is.null(data)) { 
      
      print("dimension data")
      dim(data)
      wantedData <- cleanUp(mirna1, data)
      print(ncol(wantedData))
      print(length(wantedData))
      
      if(!input$log2) wantedData[, 2] <- log2(wantedData[, 2])
      is.na(wantedData[,2]) <- sapply(wantedData[,2], is.infinite)
      
      sum_na <- sum(is.na(wantedData[,2]))
      sum_inf <- sum(is.infinite(wantedData[,2]))
      
      print("summing nA")
      print(sum_na)
      print(sum_inf)
      
      outData = list()
      outData[[1]] = wantedData
      outData[[2]] = sum_na
      return(outData)
      
    }
  })
  
  dataPredict <- reactive({
    indata <- dataCleanup();
    too_many_unknowns = 0
    outputPredict = 0
    test = NULL
    
    if (!is.null(indata)) {
      wantedData <- indata[[1]]
      sum_na <- indata[[2]]
      print(sum_na);
      
      print("unknown values:")
      print(nrow(wantedData))
      print(nrow(mirna1))
      print(0.2*nrow(mirna1))
      
      if (
        (sum_na > 0.2*nrow(mirna1))
      ) { 
        too_many_unknowns = TRUE
        outputPredict = NULL
        print("too many unknowns");
      }
      
      else {   
        
        test <- wantedData[,2]
        test <- as.numeric(test)
        
        # impute the test values
        print("imputing..")
        x <- mirna1
        new <- runImpute(cbind(x, test))
        test = new[, (ncol(x)+1) :(ncol(new)) ]
        
        # scale and normalize
        test = t(scale((test)))
        
        outputPredict <- predict(svmradial_model, test)
        print(outputPredict);
        
        # upload new data into dropbox cloud
        newData <- drop_read_csv("newData.csv")
        newData = cbind(newData[, -1], wantedData[, 2])
        newData <- data.frame(newData)
        rownames(newData) <- rownames(mirna1)
        print("new data sizes")
        print(nrow(newData))
        print(ncol(newData))
        write.csv(newData, "newData.csv", row.names = TRUE)
        drop_upload("newData.csv", dtoken=token)
        
        #retrain for every 100th entry in newData
        drop_upload("clinicalData.RData", dtoken=token)
        drop_get(path = 'clinicalData.RData', overwrite = TRUE)
      }
      
      print("svm done")
      outputResult <- list()
      
      ################## do the treatment recommendation
      
      
      
      outputResult[[1]] = outputPredict
      outputResult[[2]] = too_many_unknowns
      outputResult[[3]] = test
      
      ####### predict treatment ######################################################
      
      return(outputResult)
      
    }
    
  })
  
  ############################## treatment Recommendation ##################################
  
  treatment <- reactive({
    indata <- dataPredict();
    
    print(indata)
    
    info <- NULL
    info[3] = input$`Age`
    info[1] = input$`Gender`
    info[2] = input$`Ethnicity`
 
    treatment = NULL
    
    print("treatment 1")
    
    if (!is.null(indata)) {
      
      if (!indata[[2]]) { 
      cancer = indata[[1]]
      test = indata[[3]]
      
      ########### add the treatment recommendation function here ######################
      print("treatment 2")
      print(cancer)
      print(info)
      treatment_list = c("CESC", "ESCA", "LGG", "LUSC", "OV", "PAAD", "STAD", "TGCT", "UCS")
      if (length(which(cancer == treatment_list)) )
        treatment = recommend_three_treatments(test, as.character(cancer), c(info[[1]], info[[2]], info[[3]]))
      
      else treatment = NULL;
      }
      
    }
    
    return(treatment)
    
  })
  
  output$text1 <- renderText({ 
    
    inputResult <- dataPredict()
    treatment <- treatment()
    
    
    if (!is.null(inputResult)) { 
    
      
      
      if(inputResult[[2]] == TRUE) { 
        
        paste("Not Enough miRNA data: more than", 0.2*nrow(mirna1) , "mirna values were missing.", sep = " ");
      } else { 
        diagnosis <- inputResult[[1]]
        diagnosis_name <- cancer_name[[as.character(diagnosis)]]
        
        print("diagnosis")
        print(diagnosis)
        paste("Diagnosis: ", diagnosis, "(", diagnosis_name, ")", sep=" ")  
      }
    }
  })
  
  output$text2 <- renderUI({ 
    
    inputResult <- dataPredict()
    treatment <- treatment()
    
    if (!is.null(inputResult)) { 
      
      if(inputResult[[2]] == TRUE) { 
      } else { 
        
        if (is.null(treatment)) { 
          str1 <- paste("We are sorry. Currently, this cancer does not have a treatment recommendation in our tool")
          
        } else {
          str1 <- paste("<b>Treatment Recommendation 1:</b> ", treatment$treatments[1], " <b>with probability: </b>", round(treatment$probs[1],3)*100, sep=" ")
          strbr <- paste("")
          str2 <- paste("<b>Treatment Recommendation 2: </b> ", treatment$treatments[2], " <b>with probability: </b>", round(treatment$probs[2],3)*100, sep=" ")
          str3 <- paste("<b>Treatment Recommendation 3: </b>", treatment$treatments[3], "<b> with probability: </b>", round(treatment$probs[3],3)*100, sep=" ")
          
          HTML(paste(str1, strbr,str2, strbr, str3, strbr, sep = '<br/>'))
        }
      }
    }
    
  })
  
  
})