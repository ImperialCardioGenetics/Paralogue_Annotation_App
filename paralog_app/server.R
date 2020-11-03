library(shiny)
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)
library(writexl)

#library(tidyverse)

options(shiny.maxRequestSize=200*1024^2) #max upload size = 200 mb
enableBookmarking("url")

# options(shiny.sanitize.errors = TRUE)

shinyServer(function(input, output, session){

  get_paralog_search <- function(){
    #print(input$line)
    var<-unlist(strsplit(input$line,split="\\, |\\,|\\n|\\s|\\t"))
    input_line<-data.frame(mutation=var, stringsAsFactors = FALSE)
    input_line$mutation = stringr::str_replace_all(input_line$mutation,"[[:punct:][:space:]]","-")
    input_line$mutation = stringr::str_replace_all(input_line$mutation,"^chr","")
    input_line$paraloc = substr(input_line$mutation, 1, nchar(input_line$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
    #print(input_data)
    
    #input_line <- tidyr::separate(input_line,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
    
    #new tabix func
    result<-predict_output_tabix(validate_input(input_line))
    return(result)

  }
  
  
  get_paralog<-function(){
    # print(input$line)
    # print(input$var)
    # print(input$file)
    
    # if (input$line != "" ) {
    # 
    #   #print(input$line)
    #   var<-unlist(strsplit(input$line,split="\\, |\\,|\\n|\\s|\\t"))
    #   input_line<-data.frame(mutation=var, stringsAsFactors = FALSE)
    #   input_line$mutation = stringr::str_replace_all(input_line$mutation,"[[:punct:][:space:]]","-")
    #   input_line$mutation = stringr::str_replace_all(input_line$mutation,"^chr","")
    #   input_line$paraloc = substr(input_line$mutation, 1, nchar(input_line$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
    #   #print(input_line)
    # 
    #   #new tabix func
    #   result<-predict_output_tabix( input_line)
    #   # return(result)
    # 
    # } else 
    
    if(input$format=='paste' ){
      #shinyjs::reset("line") 
      
        #print(input$var)
        var<-unlist(strsplit(input$var,split="\\, |\\,|\\n|\\s|\\t"))
        #var=var[nzchar(x=var)]
        input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
        input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
        input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
        input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
        #print(input_data)
        
        #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
        
        #new tabix func
        result<-predict_output_tabix(validate_input(input_data))

    }else 

      if(input$format == 'upload' ) {
        #shinyjs::reset("line") 
        

        input_file <- check_upload_file(input$file)
        
        colnames(input_file) <- "mutation"
        input_file$mutation = stringr::str_replace_all(input_file$mutation,"[[:punct:][:space:]]","-")
        input_file$mutation = stringr::str_replace_all(input_file$mutation,"^chr","")
        input_file$paraloc = substr(input_file$mutation, 1, nchar(input_file$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
        
        #input_file <- tidyr::separate(input_file,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
        
        
        #new tabix func
        result<-predict_output_tabix(validate_input(input_file))

    }
        #print(input$format)
    return(result)
    
  }
    
  


  
  
  
  observeEvent((input$search_button ), {
    
    if (nrow(get_paralog_search()$paralog)>=1){ # check if result table is empty
    #if (is.na(get_paralog_search()$paralog)){ # check if result table is empty
        
      output$paralog<-renderDataTable(DT::datatable(isolate(add_paralog_URL(get_paralog_search()$paralog)),
                                                    escape = F,
                                                    extensions = c('Buttons','RowGroup'),
                                                    rownames = FALSE,
                                                    colnames = paralog_DT_colnames,
                                                    class = "display",
                                                    selection =  "none",
                                                    options = list(
                                                      rowGroup = list(dataSrc = c(5)),
                                                      #dom = 'lfrti',
                                                      dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" +
                                                             "<"row"<"col-sm-12"tr>>" +
                                                             "<"row"<"col-sm-6"i>>" +
                                                             "<"row"<"col-sm-6 btn-md"B>>"',
                                                      buttons = list(list(extend = 'excel',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                                                                          filename = paste0("paralogue_annotation_",Sys.Date())),
                                                                     list(extend = 'csv',
                                                                          fieldBoundary = '',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                                                                          fieldSeparator = '\t',
                                                                          filename = paste0("paralogue_annotation_",Sys.Date()),
                                                                          extension = '.txt')),
                                                      paging = T,scrollX = TRUE,
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(1:4,6:21,26, 30:32)),
                                                        list(width = "100px",targets = 5),
                                                        list(orderable = FALSE, className = 'details-control', targets = 0))),
                                                    callback = JS(childrow_JS_callback)
                                                    )
                                      %>% formatStyle(c(" ", "Query variant"), backgroundColor = '#f0f0f0')
                                      )

      
      output$paraloc<-renderDataTable(DT::datatable(isolate( add_paraloc_URL_new(get_paralog_search()$paraloc) ),
                                                    escape = F,
                                                    extensions = c('Buttons','RowGroup'),
                                                    rownames = FALSE,
                                                    colnames = paraloc_DT_colnames, #changed to new
                                                    class = "display",
                                                    selection =  "none",
                                                    options = list(
                                                      rowGroup = list(dataSrc = c(3)),
                                                      #dom = 'lfrti',
                                                      dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" + 
                                                             "<"row"<"col-sm-12"tr>>" + 
                                                             "<"row"<"col-sm-6"i>>" +
                                                             "<"row"<"col-sm-6 btn-md"B>>"',
                                                      buttons = list(list(extend = 'excel',
                                                                          text = ' <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                                                                          filename = paste0("paralogue_positions_",Sys.Date())),
                                                                     list(extend = 'csv',
                                                                          fieldBoundary = '',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                                                                          fieldSeparator = '\t',
                                                                          filename = paste0("paralogue_positions_",Sys.Date()),
                                                                          extension = '.txt')),
                                                      paging = T,scrollX = FALSE,
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(0:2,8,9)),
                                                        list(width = "130px",targets = 3),
                                                        list(className = 'dt-center', targets = c(3))
                                                        )
                                                      )
                                                    ) %>% 
                                        formatStyle( "Query variant", backgroundColor = '#f0f0f0')
                                      )
      
      # draw protein ----
      output$draw_prot <- renderUI({

          # fluidRow(isolate(draw_prot_data(get_paralog_search()$paralog)))
        
          fluidRow(isolate(draw_prot_data_plotly(get_paralog_search()$paralog)))
        #fluidRow(isolate(draw_prot_data_plotly(get_paralog_search()$paralog)))
      })
      
    } else {
      
      
      if (nrow(isolate(get_paralog_search()$paraloc))==0) {
        #Error catching for if query returns empty table
        isolate(showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no variants<br>Please try another input variant(s)<br>"),
          easyClose = TRUE))
        )
        
        # show NULL tables
        output$paralog<-isolate(NULL)
        output$paraloc<-isolate(NULL)
        output$draw_prot<-isolate(NULL)
        
        
        shinyjs::reset("tab2_search") 
        shinyjs::reset("All_results")
        
      } else {
        #Error catching for if query returns empty table
        updateTabsetPanel(session, "All_results", selected = "tab2")
        
        isolate(showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no paralogue variants<br>Your query has returned paralogous positions"), 
          easyClose = TRUE)))
      
        output$paralog<-isolate(NULL)
        output$draw_prot<-isolate(NULL)
        
        output$paraloc<-renderDataTable(DT::datatable(isolate( add_paraloc_URL_new(get_paralog_search()$paraloc) ),
                                                      escape = F,
                                                      extensions = c('Buttons','RowGroup'),
                                                      rownames = FALSE,
                                                      colnames = paraloc_DT_colnames, #changed to new,
                                                      class = "display",
                                                      selection =  "none",
                                                      options = list(
                                                        rowGroup = list(dataSrc = c(3)),
                                                        #dom = 'lfrti',
                                                        dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" + 
                                                             "<"row"<"col-sm-12"tr>>" + 
                                                             "<"row"<"col-sm-6"i>>" +
                                                             "<"row"<"col-sm-6 btn-md"B>>"',
                                                        buttons = list(list(extend = 'excel',
                                                                            text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                                                                            filename = paste0("paralogue_positions_",Sys.Date())),
                                                                       list(extend = 'csv',
                                                                            fieldBoundary = '',
                                                                            text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                                                                            fieldSeparator = '\t',
                                                                            filename = paste0("paralogue_positions_",Sys.Date()),
                                                                            extension = '.txt')),
                                                        paging = T,scrollX = FALSE,
                                                        columnDefs = list(
                                                          list(visible = FALSE, targets = c(0:2,8,9)),
                                                          list(width = "130px",targets = 3),
                                                          list(className = 'dt-center', targets = c(3))
                                                          )
                                                        )
                                                      ) %>% 
                                          formatStyle("Query variant", backgroundColor = '#f0f0f0')
                                        )
        }
    }
  })
    
    
    

    
  observeEvent((input$submit_button), {
    
    if (nrow(get_paralog()$paralog)>=1){ # check if result table is empty
    # if (is.na(get_paralog()$paralog)){ # check if result table is empty
        
      
      output$paralog<-renderDataTable(DT::datatable(isolate(add_paralog_URL(get_paralog()$paralog)),
                                                    escape = F,
                                                    extensions = c('Buttons','RowGroup'),
                                                    rownames = FALSE,
                                                    colnames = paralog_DT_colnames,
                                                    class = "display",
                                                    selection =  "none",
                                                    options = list(
                                                      rowGroup = list(dataSrc = c(5)),
                                                      #dom = 'lfrti',
                                                      dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" + 
                                                             "<"row"<"col-sm-12"tr>>" + 
                                                             "<"row"<"col-sm-6"i>>" +
                                                             "<"row"<"col-sm-6 btn-md"B>>"',
                                                      buttons = list(list(extend = 'excel',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                                                                          filename = paste0("paralogue_annotation_",Sys.Date())),
                                                                     list(extend = 'csv',
                                                                          fieldBoundary = '',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                                                                          fieldSeparator = '\t',
                                                                          filename = paste0("paralogue_annotation_",Sys.Date()),
                                                                          extension = '.txt')),
                                                      paging = T,scrollX = FALSE,
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(1:4,6:21,26, 30:32)),
                                                        list(width = "100px",targets = 5),
                                                        list(orderable = FALSE, className = 'details-control', targets = 0))),
                                                    callback = JS(childrow_JS_callback)
                                                    ) %>%
                                        formatStyle(c(" ", "Query variant"), backgroundColor = '#f0f0f0')
                                      )


      output$paraloc<-renderDataTable(DT::datatable(isolate( add_paraloc_URL_new(get_paralog()$paraloc) ),
                                                    escape = F,
                                                    extensions = c('Buttons','RowGroup'),
                                                    rownames = FALSE,
                                                    colnames = paraloc_DT_colnames, #changed to new
                                                    class = "display",
                                                    selection =  "none",
                                                    options = list(
                                                      rowGroup = list(dataSrc = c(3)),
                                                      #dom = 'lfrti',
                                                      dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" + 
                                                             "<"row"<"col-sm-12"tr>>" + 
                                                             "<"row"<"col-sm-6"i>>" +
                                                             "<"row"<"col-sm-6 btn-md"B>>"',
                                                      buttons = list(list(extend = 'excel',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                                                                          filename = paste0("paralogue_positions_",Sys.Date())),
                                                                     list(extend = 'csv',
                                                                          fieldBoundary = '',
                                                                          text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                                                                          fieldSeparator = '\t',
                                                                          filename = paste0("paralogue_positions_",Sys.Date()),
                                                                          extension = '.txt')),
                                                      paging = T,scrollX = FALSE,
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(0:2,8,9)),
                                                        list(width = "130px",targets = 3),
                                                        list(className = 'dt-center', targets = c(3))
                                                        )
                                                      )
                                                    ) %>%
                                        formatStyle( "Query variant", backgroundColor = '#f0f0f0')
                                      )
      

      
      # draw protein ----
      output$draw_prot <- renderUI({
        
          # fluidRow(isolate(draw_prot_data(get_paralog_search()$paralog)))
          
          fluidRow(isolate(draw_prot_data_plotly(get_paralog()$paralog)))
          
        })
      
      
    } else {
      if (nrow(isolate(get_paralog()$paraloc))==0) {
        #Error catching for if query returns empty table
        isolate(showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no variants<br>Please try another input variant(s)<br>"),
          easyClose = TRUE))
        )
        
        # show NULL tables
        output$paralog<-isolate(NULL)
        output$paraloc<-isolate(NULL)
        output$draw_prot<-isolate(NULL)
        
        shinyjs::reset("tab2_search") 
        shinyjs::reset("All_results")
        
      } else {
        #Error catching for if query returns empty table
        updateTabsetPanel(session, "All_results", selected = "tab2")
        
        isolate(showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no paralogue variants<br>Your query has returned paralogous positions"), 
          easyClose = TRUE)))
        
        output$paralog<-isolate(NULL)
        output$draw_prot<-isolate(NULL)
        

       
       output$paraloc<-renderDataTable(DT::datatable(isolate( add_paraloc_URL_new(get_paralog()$paraloc) ),
                                                     escape = F,
                                                     extensions = c('Buttons','RowGroup'),
                                                     rownames = FALSE,
                                                     colnames = paraloc_DT_colnames, #changed to new
                                                     class = "display",
                                                     selection =  "none",
                                                     options = list(
                                                       rowGroup = list(dataSrc = c(3)),
                                                       #dom = 'lfrti',
                                                       dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" + 
                                                             "<"row"<"col-sm-12"tr>>" + 
                                                             "<"row"<"col-sm-6"i>>" +
                                                             "<"row"<"col-sm-6 btn btn-md"B>>"',
                                                       buttons = list(list(extend = 'excel',
                                                                           text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                                                                           filename = paste0("paralogue_positions_",Sys.Date())),
                                                                      list(extend = 'csv',
                                                                           fieldBoundary = '',
                                                                           text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                                                                           fieldSeparator = '\t',
                                                                           filename = paste0("paralogue_positions_",Sys.Date()),
                                                                           extension = '.txt')),
                                                       paging = T,scrollX = FALSE,
                                                       columnDefs = list(
                                                         list(visible = FALSE, targets = c(0:2,8,9)),
                                                         list(width = "85px",targets = 3),
                                                         list(className = 'dt-center', targets = c(3:4))
                                                         )
                                                       )
                                                     ) %>%
                                         formatStyle("Query variant", backgroundColor = '#f0f0f0')
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
  
  observe({
    if (is.null(input$line) || input$line == "") {
    shinyjs::disable("search_button")
  } else {
    shinyjs::enable("search_button")
  }
    
  })
  
  observeEvent(input$search_button, {
    updateTabsetPanel(session, "navbar", selected = "tab2")
    #shinyjs::reset("line") 
    
  })

  observeEvent(input$submit_button, {
    #updateTabsetPanel(session, "navbar", selected = "tab2")
    shinyjs::reset("line") 
    
  })

  
  observeEvent(input$reset_button, {
    shinyjs::reset("tab2_search")
    
    
  })

  # homepage click link to tab2
  observeEvent(input$link_to_tabpanel_b, {
    updateTabsetPanel(session, "navbar", selected = "tab2")
  })
  
  
  
  
  
# Download handlers   
    
  output$download_paralog <- downloadHandler(
    filename = function() {
      paste0("paralogue_annotation_", Sys.Date(),".txt") 
    },
    content = function(file) {
      write.table(get_paralog()$paralog, file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  

  output$download_paralog_excel <- downloadHandler(
    filename = function() {
      paste0("paralogue_annotation_", Sys.Date(),".xlsx") 
    },
    content = function(file) {
      write_xlsx(x = get_paralog()$paralog, file)
      
    }
  )
  
  output$download_paraloc <- downloadHandler(
    filename = function() {
      paste0("paralogue_positions_", Sys.Date(), ".txt") 
    },
    content = function(file) {
      write.table(get_paralog()$paraloc, file, row.names = FALSE,quote = F,sep="\t")
      
    }
  )
  
  output$download_paraloc_excel <- downloadHandler(
    filename = function() {
      paste0("paralogue_positions_", Sys.Date(),".xlsx") 
    },
    content = function(file) {
      write_xlsx(x = get_paralog()$paraloc, file)
      
    }
  )
  
  
  
})






