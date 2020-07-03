library(shiny)
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)
library(writexl)

#library(tidyverse)

options(shiny.maxRequestSize=200*1024^2) #max upload size = 200 mb
enableBookmarking("url")
shinyServer(function(input, output, session){

  get_paralog_search <- function(savefile="NO"){
    # validate(
    #   need(grepl("^[[:digit:]]",input$line), "Input a valid variant"),
    #   need(nchar(input$line)<=10 , "Input very short")
    # )
    print(input$line)
    input_data<-data.frame(mutation=input$line, stringsAsFactors = FALSE)
    input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
    input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
    input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
    
    #new tabix func
    result<-predict_output_tabix(input_data)$output
    result_paraloc<-predict_output_tabix(input_data)$paraloc_output
    
    
    
    if (savefile=="NO"){
      
      if (nrow(result)!=0){ 
        result <- add_paralog_URL(result)
        result_paraloc <- add_paraloc_URL(result_paraloc)
        
      } else {
        result_paraloc <- add_paraloc_URL(result_paraloc)
      }
      
      return(list("result" = result, "result_paraloc" = result_paraloc))
      
    } else {
      
      return(list("result" = result, "result_paraloc" = result_paraloc))
    }
  }
  
  
  get_paralog<-function(savefile="NO"){
    
    if(input$format=='paste'){
      
        var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
        var=var[nzchar(x=var)]
        input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
        input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
        input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
        input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
  
        #new tabix func
        result<-predict_output_tabix(input_data)$output
        result_paraloc<-predict_output_tabix(input_data)$paraloc_output

    }else 

      if(input$format == 'upload') {

        input_file <- check_upload_file(input$file)
        
        colnames(input_file) <- "mutation"
        input_file$mutation = stringr::str_replace_all(input_file$mutation,"[[:punct:][:space:]]","-")
        input_file$mutation = stringr::str_replace_all(input_file$mutation,"^chr","")
        input_file$paraloc = substr(input_file$mutation, 1, nchar(input_file$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
        
        #new tabix func
        result<-predict_output_tabix(input_file)$output
        result_paraloc<-predict_output_tabix(input_file)$paraloc_output
      }

    if (savefile=="NO"){
      
      if (nrow(result)!=0){ 
        result <- add_paralog_URL(result)
        result_paraloc <- add_paraloc_URL(result_paraloc)

      } else {
        result_paraloc <- add_paraloc_URL(result_paraloc)
      }

      return(list("result" = result, "result_paraloc" = result_paraloc))
      
    } else {

      return(list("result" = result, "result_paraloc" = result_paraloc))
    }
  }

  
  #query_one <- reactive({get_paralog()$result})
  
  #query_positions <- reactive({get_paralog()$result_paraloc})
  
  
  observeEvent((input$search_button), {
    updateTabsetPanel(session, "navbar", selected = "tab2")
    
    if (nrow(get_paralog_search()$result)>=1){ # check if result table is empty
      #output$paralog<-renderDataTable(DT::datatable( isolate( format_DT_paralog(get_paralog()$result()) ) ) )
      
      output$paralog<-renderDataTable(DT::datatable(isolate(cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', get_paralog_search()$result)),
                                                    escape = F,
                                                    extensions = 'Buttons',
                                                    rownames = FALSE,
                                                    colnames = paralog_DT_colnames,
                                                    class = "display",
                                                    selection =  "none",
                                                    options = list(
                                                      dom = 'lfrti',
                                                      paging = T,
                                                      scrollX = TRUE,
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(1:4,6:20, 28:30)),
                                                        list(orderable = FALSE, className = 'details-control', targets = 0))
                                                      ),
                                                    callback = JS(childrow_JS_callback)
                                                    )
                                      )
      
      #output$paraloc<-renderDataTable(DT::datatable(isolate(add_paraloc_URL(get_paralog()$result_paraloc[,-c(1:3)])),
      
      output$paraloc<-renderDataTable(DT::datatable(isolate(get_paralog_search()$result_paraloc[,-c(1:3)]),
                                                    escape = F,
                                                    extensions = 'Buttons',
                                                    selection =  "none",
                                                    options = list(
                                                      dom = 'lfrti',
                                                      searchHighlight = TRUE,
                                                      paging = T, 
                                                      scrollX = FALSE,
                                                      columnDefs = list(list(width = "100px",targets = c(0)))),
                                                    rownames = FALSE,
                                                    colnames = paraloc_DT_colnames,
                                                    class = "display")
                                      )
    } else {
      
      
      if (nrow(get_paralog_search()$result_paraloc)==0) {
        
        #Error catching for if query returns empty table
        output$paralog<-showModal(modalDialog(
          title = "PARALOG Annotator", 
          HTML("Your query returned no variants<br>Please try another input variant(s)<br>"),
          easyClose = TRUE))
        
        shinyjs::reset("tab2_search") 
        
      } else {
        #Error catching for if query returns empty table
        output$paralog<-showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no paralogue variants<br>Your query has returned paralogous positions"), 
          easyClose = TRUE))
        
        updateTabsetPanel(session, "All_results", selected = "tab2")
        
        
        output$paraloc<-renderDataTable(DT::datatable(isolate(get_paralog_search()$result_paraloc[,-c(1:3)]),
                                                      #output$paraloc<-renderDataTable(DT::datatable(isolate(add_paraloc_URL(get_paralog()$result_paraloc[,-c(1:3)])),
                                                      escape = F,
                                                      extensions = 'Buttons',
                                                      options = list(
                                                        dom = 'lfrti',
                                                        searchHighlight = TRUE,
                                                        paging = T, 
                                                        scrollX = FALSE,
                                                        columnDefs = list(list(width = "100px",targets = c(0)))),
                                                      rownames = FALSE,
                                                      colnames = c(
                                                        'Query variant' = 'var.query',
                                                        'Gene' = 'Gene.query',
                                                        'ENSG' = 'ENSG.query',
                                                        'Paralogous positions' = 'Positions.paralog'),
                                                      class = "display")
                                        )
        }
    }
    #shinyjs::reset("All_results")
  })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  observeEvent((input$submit_button ), {
    
    if (nrow(get_paralog()$result)>=1){ # check if result table is empty
      #output$paralog<-renderDataTable(DT::datatable( isolate( format_DT_paralog(get_paralog()$result()) ) ) )
      
      output$paralog<-renderDataTable(DT::datatable(isolate(cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', get_paralog()$result)),
                                                    escape = F,
                                                    extensions = 'Buttons',
                                                    rownames = FALSE,
                                                    colnames = paralog_DT_colnames,
                                                    class = "display",
                                                    selection =  "none",
                                                    options = list(
                                                      dom = 'lfrti',
                                                      paging = T,
                                                      scrollX = TRUE,
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(1:4,6:20, 28:30)),
                                                        list(orderable = FALSE, className = 'details-control', targets = 0))
                                                      ),
                                                    callback = JS(childrow_JS_callback)
                                                    )
                                      )

      #output$paraloc<-renderDataTable(DT::datatable(isolate(add_paraloc_URL(get_paralog()$result_paraloc[,-c(1:3)])),
                                                    
      output$paraloc<-renderDataTable(DT::datatable(isolate(get_paralog()$result_paraloc[,-c(1:3)]),
                                                    escape = F,
                                                    extensions = 'Buttons',
                                                    selection =  "none",
                                                    options = list(
                                                      dom = 'lfrti',
                                                      searchHighlight = TRUE,
                                                      paging = T, 
                                                      scrollX = FALSE,
                                                      columnDefs = list(list(width = "100px",targets = c(0)))),
                                                    rownames = FALSE,
                                                    colnames = paraloc_DT_colnames,
                                                    class = "display")
                                      )
    } else {
      
      
      if (nrow(get_paralog()$result_paraloc)==0) {
        
        #Error catching for if query returns empty table
        output$paralog<-showModal(modalDialog(
          title = "PARALOG Annotator", 
          HTML("Your query returned no variants<br>Please try another input variant(s)<br>"),
          easyClose = TRUE))
        
        shinyjs::reset("tab2_search") 
      
      } else {
        #Error catching for if query returns empty table
        output$paralog<-showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no paralogue variants<br>Your query has returned paralogous positions"), 
          easyClose = TRUE))
        
       updateTabsetPanel(session, "All_results", selected = "tab2")
       
       
       output$paraloc<-renderDataTable(DT::datatable(isolate(get_paralog()$result_paraloc[,-c(1:3)]),
       #output$paraloc<-renderDataTable(DT::datatable(isolate(add_paraloc_URL(get_paralog()$result_paraloc[,-c(1:3)])),
                                                     escape = F,
                                                     extensions = 'Buttons',
                                                     options = list(
                                                       dom = 'lfrti',
                                                       searchHighlight = TRUE,
                                                       paging = T, 
                                                       scrollX = FALSE,
                                                       columnDefs = list(list(width = "100px",targets = c(0)))),
                                                     rownames = FALSE,
                                                     colnames = c(
                                                       'Query variant' = 'var.query',
                                                       'Gene' = 'Gene.query',
                                                       'ENSG' = 'ENSG.query',
                                                       'Paralogous positions' = 'Positions.paralog'),
                                                     class = "display")
                                       )
      }
    }
  })
  
  observe({
    if (is.null(input$var) || input$var == ""){
      if (is.null(input$file) || input$file == ""){
        shinyjs::disable("submit_button")
        shinyjs::disable("reset_button")
      } else {
        shinyjs::enable("submit_button")
        shinyjs::enable("reset_button")
      }
    } else {
      shinyjs::enable("submit_button")
      shinyjs::enable("reset_button")
    }
  })
  
  observe({if (is.null(input$line) || input$line == "") {
    shinyjs::disable("search_button")
  } else {
    shinyjs::enable("search_button")
  }
    
  })
  
  # observeEvent(input$search_button, {
  #   updateTabsetPanel(session, "navbar", selected = "tab2")
  # })
  # 
  # 
  
  observeEvent(input$reset_button, {
    shinyjs::reset("tab2_search")
    
    
  })

  
  
  
  
  
  
# Download handlers   
    
  output$download_paralog <- downloadHandler(
    filename = function() {
      paste0("paralogue_annotation_", Sys.Date(),".tsv") 
    },
    content = function(file) {
      write.table(get_paralog(savefile = "YES")$result, file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  

  output$download_paralog_excel <- downloadHandler(
    filename = function() {
      paste0("paralogue_annotation_", Sys.Date(),".xlsx") 
    },
    content = function(file) {
      write_xlsx(x = get_paralog(savefile = "YES")$result, file)
      
    }
  )
  
  output$download_paraloc <- downloadHandler(
    filename = function() {
      paste0("paralogue_positions_", Sys.Date(), ".tsv") 
    },
    content = function(file) {
      write.table(get_paralog(savefile = "YES")$result_paraloc, file, row.names = FALSE,quote = F,sep="\t")
      
    }
  )
  
  output$download_paraloc_excel <- downloadHandler(
    filename = function() {
      paste0("paralogue_positions_", Sys.Date(),".xlsx") 
    },
    content = function(file) {
      write_xlsx(x = get_paralog(savefile = "YES")$result_paraloc, file)
      
    }
  )
  
  
})






