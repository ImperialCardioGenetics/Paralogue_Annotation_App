suppressWarnings(suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinyjs)
  library(shinycssloaders)
  library(tidyverse)
  library(DT)
  library(kableExtra)
}))

options(shiny.maxRequestSize=10*1024^2) #max upload size = 100 mb
enableBookmarking("url")

# options(shiny.sanitize.errors = TRUE)

shinyServer(function(input, output, session){
  
  # generate Paraloc RAW data Download Table
  output$table <- renderUI({generate_download_table(output, session)
  })
  
  
  # set search and main functionss for tabs separately
  get_paralog_search <- function(){lookup_paralog_new(validate_input(input$line))}
  get_paraloc_search <- function(){lookup_paraloc_new(validate_input(input$line))}
  get_homolog_search <- function(){lookup_homolog_new(validate_input(input$line))}

  get_paralog_main <- function(){if (input$format=='paste') {lookup_paralog_new(validate_input(input$var))} else {lookup_paralog_new(validate_input(check_upload_file(input$file)))}}
  get_paraloc_main <- function(){if (input$format=='paste') {lookup_paraloc_new(validate_input(input$var))} else {lookup_paraloc_new(validate_input(check_upload_file(input$file)))}}
  get_homolog_main <- function(){if (input$format=='paste') {lookup_homolog_new(validate_input(input$var))} else {lookup_homolog_new(validate_input(check_upload_file(input$file)))}}
    
  
  # get number of vars output for modal msg
  n_get_paralog_search <- reactive({nrow(get_paralog_search())})
  n_get_paraloc_search <- reactive({nrow(get_paraloc_search())})
  n_get_homolog_search <- reactive({nrow(get_homolog_search())})
  
  n_get_paralog_main   <- reactive({nrow(get_paralog_main())})
  n_get_paraloc_main   <- reactive({nrow(get_paraloc_main())})
  n_get_homolog_main   <- reactive({nrow(get_homolog_main())})
  

  # write DT databale object
  paralog_search<- reactive({renderDataTable(DT::datatable(isolate(add_URL(get_paralog_search())),
                                                           escape = F,
                                                           extensions = c('Buttons','RowGroup'),
                                                           rownames = FALSE,
                                                           colnames = paralog_DT_colnames,
                                                           class = "display",
                                                           selection =  "none",
                                                           options = paralog_DT_options_list,
                                                           callback = JS(childrow_JS_callback_paralog)) %>% formatStyle(c(" ", "Query variant"),backgroundColor = '#f0f0f0'))
    })
  
  paralog<- reactive({renderDataTable(DT::datatable(isolate(add_URL(get_paralog_main())), 
                                                    escape = F,
                                                    extensions = c('Buttons','RowGroup'),
                                                    rownames = FALSE,
                                                    colnames = paralog_DT_colnames,
                                                    class = "display",
                                                    selection =  "none",
                                                    options = paralog_DT_options_list,
                                                    callback = JS(childrow_JS_callback_paralog)) %>% formatStyle(c(" ", "Query variant"),backgroundColor = '#f0f0f0'))

    })
  

  paraloc_search <- reactive({renderDataTable(DT::datatable(isolate(add_URL(get_paraloc_search())),
                                                            escape = F,
                                                            extensions = c('Buttons','RowGroup'),
                                                            rownames = FALSE,
                                                            colnames = paraloc_DT_colnames, #changed to new
                                                            class = "display",
                                                            selection =  "none",
                                                            options = paraloc_DT_options_list) %>% formatStyle( c("Query variant","Query gene", "Query residue"), backgroundColor = '#f0f0f0'))
    })
  

  paraloc <- reactive({renderDataTable(DT::datatable(isolate(add_URL(get_paraloc_main())),
                                                     escape = F,
                                                     extensions = c('Buttons','RowGroup'),
                                                     rownames = FALSE,
                                                     colnames = paraloc_DT_colnames, #changed to new
                                                     class = "display",
                                                     selection =  "none",
                                                     options = paraloc_DT_options_list) %>% formatStyle( c("Query variant","Query gene", "Query residue"), backgroundColor = '#f0f0f0'))
    
  })
  
  
  homolog_search<- reactive({renderDataTable(DT::datatable(isolate(add_URL(get_homolog_search())),
                                                           escape = F,
                                                           extensions = c('Buttons','RowGroup'),
                                                           rownames = FALSE,
                                                           colnames = homolog_DT_colnames,
                                                           class = "display",
                                                           selection =  "none",
                                                           options = homolog_DT_options_list,
                                                           callback = JS(childrow_JS_callback_homolog)) %>% formatStyle(c(" ", "Query variant"),backgroundColor = '#f0f0f0'))
  })
  

  
  homolog <- reactive({renderDataTable(DT::datatable(isolate(add_URL(get_homolog_main())),
                                                     escape = F,
                                                     extensions = c('Buttons','RowGroup'),
                                                     rownames = FALSE,
                                                     colnames = homolog_DT_colnames,
                                                     class = "display",
                                                     selection =  "none",
                                                     options = homolog_DT_options_list,
                                                     callback = JS(childrow_JS_callback_homolog)) %>% formatStyle(c(" ", "Query variant"),backgroundColor = '#f0f0f0'))

  })
  
  
  
  observeEvent(input$search_button,{
    

    #Error catching for if query returns empty table
    #if (nrow(get_paralog_search()$paralog)>=1) {
    if (n_get_paralog_search()>=1) {
        
      
      if (n_get_homolog_search()>=1) {
        
        output$paralog <- isolate(paralog_search())
        output$paraloc <- isolate(paraloc_search())
        output$homolog <- isolate(homolog_search())
        
        } else {
        
        # display modal
        isolate(showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no homologous pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned paralogous positions and paralogous missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b>"), 
          easyClose = TRUE)))
        
        output$paralog<- isolate(paralog_search())
        output$paraloc <- isolate(paraloc_search())
        output$homolog <- isolate(NULL)
        
        }
      } else {
        #Error catching for if query returns empty table
        if (n_get_paraloc_search()>=1) {
          
          if (n_get_homolog_search()>=1) {
            
            updateTabsetPanel(session, "All_results", selected = "tab2")
            
            # display modal
            isolate(showModal(modalDialog(
              title = "PARALOG Annotator",
              HTML("Your query returned no paralogous missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned paralogous positions and homologous pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b>"), 
              easyClose = TRUE)))
            
            output$paralog <- isolate(NULL)
            output$paraloc <- isolate(paraloc_search())
            output$homolog <- isolate(homolog_search())
            
            } else {
              
              updateTabsetPanel(session, "All_results", selected = "tab2")
              
              # display modal
              isolate(showModal(modalDialog(
                title = "PARALOG Annotator",
                HTML("Your query returned no paralogous or homologous Pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned paralogous positions"), 
                easyClose = TRUE)))
              
              output$paralog <- isolate(NULL)
              output$paraloc <- isolate(paraloc_search())
              output$homolog <- isolate(NULL)
              
              }
          
          } else {
            
            if (n_get_homolog_search()>=1) {
              
              updateTabsetPanel(session, "All_results", selected = "tab3")
              
              # display modal
              isolate(showModal(modalDialog(
                title = "PARALOG Annotator",
                HTML("Your query returned no paralogous positions or paralogous missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned homologous Pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b>"), 
                easyClose = TRUE)))
              
              output$paralog <- isolate(NULL)
              output$paraloc <- isolate(NULL)
              output$homolog <- isolate(homolog_search())
              
              } else {
              
                # display modal
                isolate(showModal(modalDialog(
                  title = "PARALOG Annotator",
                  HTML("Your query returned no results<br><br>Please try another input variant(s)<br>"),
                  easyClose = TRUE))
                )
                
                # show NULL tables
                output$paralog <- isolate(NULL)
                output$paraloc <- isolate(NULL)
                output$homolog <- isolate(NULL)
                
                
                shinyjs::reset("tab2_search") 
                shinyjs::reset("All_results")
                
                }
          
          }
      }
  })
    
    
    

    
  observeEvent(input$submit_button, {
    
    #Error catching for if query returns empty table
    if (n_get_paralog_main()>=1) {
      
      if (n_get_homolog_main()>=1) {
        
        output$paralog <- isolate(paralog())
        output$paraloc <- isolate(paraloc())
        output$homolog <- isolate(homolog())
        
      } else {
        
        # display modal
        isolate(showModal(modalDialog(
          title = "PARALOG Annotator",
          HTML("Your query returned no homologous pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned paralogous positions and paralogous missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b>"), 
          easyClose = TRUE)))
        
        output$paralog <- isolate(paralog())
        output$paraloc <- isolate(paraloc())
        output$homolog <- isolate(NULL)
        
      }
      
    } else {
      
      #Error catching for if query returns empty table
      if (n_get_paraloc_main()>=1) {
        
        if (n_get_homolog_main()>=1) {
          
          updateTabsetPanel(session, "All_results", selected = "tab2")
          
          # display modal
          isolate(showModal(modalDialog(
            title = "PARALOG Annotator",
            HTML("Your query returned no paralogous missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned paralogous positions and homologous pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b>"), 
            easyClose = TRUE)))
          
          output$paralog <- isolate(NULL)
          output$paraloc <- isolate(paraloc())
          output$homolog <- isolate(homolog())
          
        } else {
          
          updateTabsetPanel(session, "All_results", selected = "tab2")
          
          # display modal
          isolate(showModal(modalDialog(
            title = "PARALOG Annotator",
            HTML("Your query returned no paralogous or homologous Pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned paralogous positions"), 
            easyClose = TRUE)))
          
          output$paralog <- isolate(NULL)
          output$paraloc <- isolate(paraloc())
          output$homolog <- isolate(NULL)
          
        }
        
      } else {
        
        if (n_get_homolog_main()>=1) {
          
          updateTabsetPanel(session, "All_results", selected = "tab3")
          
          # display modal
          isolate(showModal(modalDialog(
            title = "PARALOG Annotator",
            HTML("Your query returned no paralogous positions or paralogous missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b><br><br>Your query has returned homologous Pfam missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b>"), 
            easyClose = TRUE)))
          
          output$paralog <- isolate(NULL)
          output$paraloc <- isolate(NULL)
          output$homolog <- isolate(homolog())
          
        } else {
          
          # display modal
          isolate(showModal(modalDialog(
            title = "PARALOG Annotator",
            HTML("Your query returned no results<br><br>Please try another input variant(s)<br>"),
            easyClose = TRUE))
          )
          
          # show NULL tables
          output$paralog <- isolate(NULL)
          output$paraloc <- isolate(NULL)
          output$homolog <- isolate(NULL)
          
          
          shinyjs::reset("tab2_search") 
          shinyjs::reset("All_results")
          
        }
        
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
    shinyjs::reset("var") 
    shinyjs::reset("file") 
    
    
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
  
  rmarkdown::render("README.md",
                    output_format = rmarkdown::html_fragment(),
                    output_file = "www/README.html",
                    output_options = list(metadata = list(title = "About")),
                    quiet = TRUE)

  # output$about<-renderUI({includeHTML("about.html")})
  output$about<-renderUI({
    includeHTML("www/README.html")
    })
  
  

  
  # download data
  output$download_file <- downloadHandler(
    filename = function() {
      "paralog_data.txt.gz"  # The name that will appear to the user
    },
    content = function(file) {
      file.copy("paralog_data.txt.gz", file)  # Replace with the actual file path
    }
  )
  
})






