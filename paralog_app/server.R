library(shiny)
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)
#library(WriteXLS)
library(writexl)
#library(tidyverse)
options(shiny.maxRequestSize=300*1024^2) #max upload size = 200 mb
enableBookmarking("url")
shinyServer(function(input, output, session){
  output$text1 <- renderText({ paste("WARNING: only Chromosome 21 available in current DEMO version",input$n) })
  # observe({
  #   # Trigger this observer every time an input changes
  #   reactiveValuesToList(input)
  #   session$doBookmark()
  # })
  # onBookmarked(function(url) {
  #   updateQueryString(url)
  # })
  
  get_paralog<-function(savefile="NO"){

      if(input$format=='paste'){
      #input<-data.frame(var="1:114713907:T:G",stringsAsFactors = F)  
      #input$var<-data.frame(var="1\t114713907\tT\tG")
        #req(input$var)
        # print(input$var)
        var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
        # print(var)
        var=var[nzchar(x=var)]
        input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
        input_data$mutation = stringr::str_replace_all(input_data$mutation,":"," ")
        input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
        input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
  
        #new tabix func
        result<-predict_output_tabix(input_data)$output
        result_paraloc<-predict_output_tabix(input_data)$paraloc_output
        
        # result<-predict_output(input_data)$output
        # result_paraloc<-predict_output(input_data)$paraloc_output
    }else{
      if(input$format == 'upload') {
        #req(input$file)
        inFile <- input$file
        
        # very hacky way to read in vcf
        # read 1st line only to get number of cols to check if its txt or vcf format
        # when using colClasses in read.table the columns set to NULL are completely ignored
        # input_1row = ncol(read.table(inFile$datapath,nrows = 1 ))
        
        input_file <- check_upload_file(inFile)
        
        
        colnames(input_file) <- "mutation"
        input_file$mutation = stringr::str_replace_all(input_file$mutation,":"," ")
        input_file$mutation = stringr::str_replace_all(input_file$mutation,"^chr","")
        input_file$paraloc = substr(input_file$mutation, 1, nchar(input_file$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
        
        #new tabix func
        result<-predict_output_tabix(input_file)$output
        result_paraloc<-predict_output_tabix(input_file)$paraloc_output
        
        
        # result <- predict_output(input_file)$output
        # result_paraloc<-predict_output(input_file)$paraloc_output
      }
    }
  #}
    if (savefile=="NO"){
      #add ClinVar IDs with URLs 
      #result$Query_ClinVar_link<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$Query_ClinVar,"/"), "' target='_blank'>", result$Query_ClinVar, "</a>")  #Not possible for custom_ids; Can add feature to check P/LP tableized file
      if (nrow(result)!=0){ # that where the error was generated
        
        #ClinVarID paralog URL
        result$ID.paralog<- ifelse(#!is.na(result$ID.paralog),
                                   result$ID.paralog!="NA",
                                   (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.paralog,"/"), "' target='_blank'>", result$ID.paralog, "</a>")),
                                   "-")
        
        #ClinVarID query URL
        result$ID.query<- ifelse(#!is.na(result$ID.query),
                                 result$ID.query!="NA",
                                 (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.query,"/"), "' target='_blank'>", result$ID.query, "</a>")),
                                 "-")
        #print(paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene.query)],";g1=",map[unlist(result$SYMBOL.paralog)]))
        #Ensembl alignment URL
        # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
        result$Ensembl_alignment_link<- ifelse(!is.na(result$SYMBOL.paralog), 
                                               (paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene.query)],";g1=",map[unlist(result$SYMBOL.paralog)]), "' class='btn btn-default btn-sm btn-block active' target='_blank'>alignment</a>")) , 
                                               "-") 
        
        # https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000213281;t=ENST00000369535
        #Ensembl Transcript.query
        result$Transcript.query<- ifelse(!is.na(result$Transcript.query), 
                                         (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",map[unlist(result$Gene.query)],";t=",result$Transcript.query), "' target='_blank'>", result$Transcript.query, "</a>")),
                                         "-")
        
        # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000213281
        #Ensembl Gene.query
        result$Gene.query<- ifelse(!is.na(result$Gene.query), 
                                   (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result$Gene.query)]), "' target='_blank'>", result$Gene.query, "</a>")),
                                   "-")
        
        #Ensembl SYMBOL.paralog
        result$SYMBOL.paralog<- ifelse(!is.na(result$SYMBOL.paralog), 
                                       (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result$SYMBOL.paralog)]), "' target='_blank'>", result$SYMBOL.paralog, "</a>")),
                                       
                                       "-")
        
        #Ensembl Gene.query for paraloc
        result_paraloc$Gene.query<- ifelse(!is.na(result_paraloc$Gene.query), 
                                   (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene.query)]), "' target='_blank'>", result_paraloc$Gene.query, "</a>")),
                                   "-")

      } else {
        #Ensembl Gene.query for paraloc
        result_paraloc$Gene.query<- ifelse(!is.na(result_paraloc$Gene.query), 
                                           (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene.query)]), "' target='_blank'>", result_paraloc$Gene.query, "</a>")),
                                           "-")
      }

      
      return(list("result" = result, "result_paraloc" = result_paraloc))
    } else {
      

      
      return(list("result" = result, "result_paraloc" = result_paraloc))
    }
  }

  
  #query_one <- reactive({get_paralog()$result})
  
  #query_positions <- reactive({get_paralog()$result_paraloc})
  
  
  observeEvent(input$submit_button, {
    
    
    
    # updateQueryString(paste0("?",input$mutation[1]), mode = "push")
    if (nrow(get_paralog()$result)>=1){ # check if result table is empty
      output$paralog<-renderDataTable(DT::datatable(isolate(cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', get_paralog()$result)),
                                                    escape = F,
                                                    extensions = 'Buttons',
                                                    rownames = FALSE,
                                                    colnames = c('Chr' = 'CHR.query',
                                                                 'Pos' = 'POS.query',
                                                                 'REF' = 'REF.query',
                                                                 'ALT' = 'ALT.query',
                                                                 ' ' = 'var.query',
                                                                 'ClinVar ID' = 'ID.query',
                                                                 'Gene' = 'Gene.query',
                                                                 'Transcript' = 'Transcript.query', 
                                                                 'Protein' = 'Protein_dot.query', 
                                                                 'Codons' = 'Codons.query',
                                                                 'Para_Z Score'='Para_Z_score.query',
                                                                 'Chr' = 'CHR.paralog',
                                                                 'Pos' = 'POS.paralog',
                                                                 'REF' = 'REF.paralog',
                                                                 'ALT' = 'ALT.paralog',
                                                                 'Paralogue variant' = 'var.paralog',
                                                                 'ClinVar ID' = 'ID.paralog',
                                                                 'Gene' = 'SYMBOL.paralog',
                                                                 'Protein' = 'Protein_dot.paralog',
                                                                 'Codon' = 'Codons.paralog',
                                                                 'Para_Z Score' = 'Para_Z_score.paralog',
                                                                 'Ensembl alignment' = 'Ensembl_alignment_link'
                                                                 ),
                                                    class = "display",
                                                    selection =  "none",
                                                    #container = sketch, 
                                                    options = list(
                                                      dom = 'lfrti',
                                                      buttons = list('copy', list(extend = 'collection',buttons = list(list(extend = 'excel',
                                                                                                                            filename = 'PARALOG_Annotator'),
                                                                                                                       list(extend = 'csv',
                                                                                                                            fieldBoundary = '',
                                                                                                                            text = 'TXT',
                                                                                                                            fieldSeparator = '\t',
                                                                                                                            filename = 'PARALOG_Annotator',
                                                                                                                            extension = '.txt'),
                                                                                                                       list(extend = 'pdf',
                                                                                                                            pageSize = 'A4',
                                                                                                                            orientation = 'landscape',
                                                                                                                            filename = 'PARALOG_Annotator')),
                                                                                  text = 'Download')),
                                                      paging = T,
                                                      scrollX = TRUE, 
                                                      columnDefs = list(
                                                        list(visible = FALSE, targets = c(1:4,6:15)), # was 2:7 with old_func
                                                        list(orderable = FALSE, className = 'details-control', targets = 0)
                                                      )
                                                    ),
                                                    callback = JS("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<table style=\"border-spacing:50px\" style=\"width:50%\" cellpadding=\"50px\" style=\"padding-left:50px\">'+
      '<thead>'+
    //  '<tr>'+
    //        '<th colspan=\"4\">Query variant</td>'+
    //  '</tr>'+
        '<tr>'+
    //      '<th>Chr Pos REF ALT</th>'+
            '<th>Query variant</th>'+
            '<th>ClinVar ID</th>'+
            '<th>HGVS</th>'+
            '<th>Codons</th>'+
            '<th>Para_Z Score</th>'+
        '</tr>'+
       '</thead>'+
       '<tbody>'+
        '<tr>'+
          //  '<td>'+ d[1] +'</td>'+
          //  '<td>'+ d[2] +'</td>'+
          //  '<td>'+ d[4] + '(' + d[3] + ').cDNA (' + d[5] + ')' +'</td>'+
          //  '<td>'+ d[6] +'</td>'+
          //  '<td>'+ d[7] +'</td>'+
            '<td>'+ d[1] +' '+ d[2] + ' '+ d[3] + ' '+ d[4] +'</td>'+
            '<td>'+ d[6] +'</td>'+
            '<td>'+ d[9] + '(' + d[7] + ').cDNA (' + d[10] + ')' +'</td>'+
            '<td>'+ d[8] +'</td>'+
            '<td>'+ d[11] +'</td>'+
        '</tr>'+
      '</tbody>'+
    '</table>';
  };
  table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
      td.html('<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css\"> <i class=\"fa fa-plus-square fa-lg\"></i>');
    } else {
      row.child(format(row.data())).show();
      td.html('<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css\"> <i class=\"fa fa-minus-square fa-lg\"></i>');
    }
  }
);"
                                                    )
                                                    )
                                      )
 

      output$paraloc<-renderDataTable(DT::datatable(isolate(add_paraloc_URL(get_paralog()$result_paraloc[,-c(1:3)])),
                                                    escape = F, # escape text hyperlink to url instead of text
                                                    extensions = 'Buttons',
                                                    options = list(
                                                      dom = 'lfrti',
                                                      buttons = list('copy', list(extend = 'collection',buttons = list(list(extend = 'excel',
                                                                                                                            filename = 'PARALOG_Annotator_locations'),
                                                                                                                       list(extend = 'csv',
                                                                                                                            fieldBoundary = '',
                                                                                                                            text = 'TXT',
                                                                                                                            fieldSeparator = '\t',
                                                                                                                            filename = 'PARALOG_Annotator_locations',
                                                                                                                            extension = '.txt'),
                                                                                                                       list(extend = 'pdf',
                                                                                                                            pageSize = 'A4',
                                                                                                                            orientation = 'landscape',
                                                                                                                            filename = 'PARALOG_Annotator_locations')),
                                                                                  text = 'Download')),
                                                      searchHighlight = TRUE,
                                                      paging = T, # True to activate pagination
                                                      scrollX = FALSE,
                                                      columnDefs = list(list(width = "100px",targets = c(0)))),# set options for table eg. per page lines #,columnDefs = list(list(className = 'dt-right', targets = c(1,5,7,11)))
                                                    rownames = FALSE,
                                                    colnames = c(#'Query position' = 'var',
                                                      #'Chr' = 'CHR.query',
                                                      #'Pos' = 'POS.query',
                                                      #'REF'= 'REF.query',
                                                      'Query variant' = 'var.query',
                                                      'Gene' = 'Gene.query',
                                                      'Paralogous positions' = 'Positions.paralog'),
                                                    class = "display")
                                      )
    } else {
      
      
      if (nrow(get_paralog()$result_paraloc)==0) {
        
        #Error catching for if query returns empty table
        output$paralog<-showModal(modalDialog(
          title = "PARALOG Annotator", # We can change the msg
          HTML("Your query returned no variants<br>Please try another input variant(s)<br>"), # and this msg
          easyClose = TRUE))
        
        shinyjs::reset("myapp") # we can delete this so the app does not restart every time
      
      } else {
        #Error catching for if query returns empty table
        output$paralog<-showModal(modalDialog(
          title = "PARALOG Annotator", # We can change the msg
          HTML("Your query returned no paralogue variants<br>Your query has returned paralogous positions"), # and this msg
          easyClose = TRUE))
        
       updateTabsetPanel(session, "All_results", selected = "tab2") # STARTED FIXING THIS , have not found a way to checnge to a conditional tab yet.
       
       output$paraloc<-renderDataTable(DT::datatable(isolate(add_paraloc_URL(get_paralog()$result_paraloc[,-c(1:3)])),
                                                     escape = F, # escape text hyperlink to url instead of text
                                                     extensions = 'Buttons',
                                                     options = list(
                                                       dom = 'lfrti',
                                                       buttons = list('copy', list(extend = 'collection',buttons = list(list(extend = 'excel',
                                                                                                                             filename = 'PARALOG_Annotator_locations'),
                                                                                                                        list(extend = 'csv',
                                                                                                                             fieldBoundary = '',
                                                                                                                             text = 'TXT',
                                                                                                                             fieldSeparator = '\t',
                                                                                                                             filename = 'PARALOG_Annotator_locations',
                                                                                                                             extension = '.txt'),
                                                                                                                        list(extend = 'pdf',
                                                                                                                             pageSize = 'A4',
                                                                                                                             orientation = 'landscape',
                                                                                                                             filename = 'PARALOG_Annotator_locations')),
                                                                                   text = 'Download')),
                                                       searchHighlight = TRUE,
                                                       paging = T, # True to activate pagination
                                                       scrollX = FALSE,
                                                       columnDefs = list(list(width = "100px",targets = c(0)))),# set options for table eg. per page lines #,columnDefs = list(list(className = 'dt-right', targets = c(1,5,7,11)))
                                                     rownames = FALSE,
                                                     colnames = c(#'Query position' = 'var',
                                                       #'Chr' = 'CHR.query',
                                                       #'Pos' = 'POS.query',
                                                       #'REF'= 'REF.query',
                                                       'Query variant' = 'var.query',
                                                       'Gene' = 'Gene.query',
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
  
  observeEvent(input$reset_button, {
    shinyjs::reset("myapp")
    # output$paralog<-renderText(isolate({
    #   
    # }))
  })
  
  output$download_paralog <- downloadHandler(
    filename = function() {
      paste0("paralogue_annotation_", Sys.Date(),".tsv") # need to give specific name?
    },
    content = function(file) {
      write.table(edit_download_cols(get_paralog(savefile = "YES")$result), file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  
  output$download_paralog_excel <- downloadHandler(
    filename = function() {
      paste0("paralogue_annotation_", Sys.Date(),".xlsx") # need to give specific name?
    },
    content = function(file) {
      #WriteXLS(x = (edit_download_cols(get_paralog(savefile = "YES")$result)), ExcelFileName = file, SheetNames = "1" )
      write_xlsx(x = edit_download_cols(get_paralog(savefile = "YES")$result), file)
      
    }
  )
  
  output$download_paraloc <- downloadHandler(
    filename = function() {
      paste0("paralogue_positions_", Sys.Date(), ".tsv") # need to give specific name?
    },
    content = function(file) {
      #write.table(edit_download_cols_paraloc(get_paralog(savefile = "YES")$result_paraloc), file, row.names = FALSE,quote = F,sep="\t")
      write.table(get_paralog(savefile = "YES")$result_paraloc, file, row.names = FALSE,quote = F,sep="\t")
      
    }
  )
  
  output$download_paraloc_excel <- downloadHandler(
    filename = function() {
      paste0("paralogue_positions_", Sys.Date(),".xlsx") # need to give specific name?
    },
    content = function(file) {
      #WriteXLS(x = (get_paralog(savefile = "YES")$result_paraloc), ExcelFileName = file, SheetNames = "1")
      write_xlsx(x = get_paralog(savefile = "YES")$result_paraloc, file)
      
    }
  )
  
  
})






