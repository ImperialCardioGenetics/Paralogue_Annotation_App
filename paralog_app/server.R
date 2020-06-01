library(shiny)
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)

#library(tidyverse)
options(shiny.maxRequestSize=200*1024^2) #max upload size = 200 mb
enableBookmarking("url")
shinyServer(function(input, output, session){
  
  # observe({
  #   # Trigger this observer every time an input changes
  #   reactiveValuesToList(input)
  #   session$doBookmark()
  # })
  # onBookmarked(function(url) {
  #   updateQueryString(url)
  # })
  
  get_paralog<-function(savefile){
    
    #input<-data.frame(chr="1",pos="114713907",ref="T",alt="A")
    # if(input$format=='pick'){
    #   req(input$chr)
    #   req(input$pos)
    #   req(input$ref)
    #   req(input$alt)
    #   var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
    #   var = var[nzchar(x=var)]
    #   input_data = data.frame(mutation=var, stringsAsFactors = FALSE)
    #   # input_data = data.frame(chr = input$chr, pos = input$pos, ref = input$ref, alt = input$alt) #not needed
    #   result<-predict_output(input_data)
    # }else{
      if(input$format=='paste'){
      #input<-data.frame(var="1:114713907:T:G",stringsAsFactors = F)  
      #input$var<-data.frame(var="1\t114713907\tT\tG")
        req(input$var)
        # print(input$var)
        var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
        # print(var)
        var=var[nzchar(x=var)]
        input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
        input_data$mutation = stringr::str_replace_all(input_data$mutation,":"," ")
        input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
        input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2) #CAN OPTIMISED THIS MAYBE LATER, JUST GETTING IT TO WORK FOR NOW
        # print(input_data)
        # colnames(input_data)<-"mutation"
        result<-predict_output(input_data)$output
        result_paraloc<-predict_output(input_data)$paraloc_output
    }else{
      if(input$format == 'upload') {
        req(input$file)
        inFile <- input$file
        
        # very hacky way to read in vcf
        # read 1st line only to get number of cols to check if its txt or vcf format
        # when using colClasses in read.table the columns set to NULL are completely ignored
        # input_1row = ncol(read.table(inFile$datapath,nrows = 1 ))
        
        input_file <- check_upload_file(inFile)
        
        
        colnames(input_file) <- "mutation"
        input_file$mutation<-gsub(":|\t"," ",input_file$mutation)
        input_file$mutation<-gsub("^chr","",input_file$mutation)
        result <- predict_output(input_file)$output
        result_paraloc<-predict_output(input_file)$paraloc_output
      }
    }
  #}
    if (savefile=="NO"){
      #add ClinVar IDs with URLs 
      #result$Query_ClinVar_link<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$Query_ClinVar,"/"), "' target='_blank'>", result$Query_ClinVar, "</a>")  #Not possible for custom_ids; Can add feature to check P/LP tableized file
      if (nrow(result)!=0){ # that where the error was generated
        
        #ClinVarID paralog URL
        result$ID.paralog<- ifelse(!is.na(result$ID.paralog), 
                                   (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.paralog,"/"), "' target='_blank'>", result$ID.paralog, "</a>")),
                                   "-")
        
        #ClinVarID query URL
        result$ID.query<- ifelse(!is.na(result$ID.query), 
                                 (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.query,"/"), "' target='_blank'>", result$ID.query, "</a>")),
                                 "-")
        #print(paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene.query)],";g1=",map[unlist(result$SYMBOL.paralog)]))
        #Ensembl alignment URL
        # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
        result$Ensembl_alignment_link<- ifelse(!is.na(result$SYMBOL), 
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

        #Ensembl Paraloc_data$Gene
        result_paraloc$Gene<- ifelse(!is.na(result_paraloc$Gene), 
                                   (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene)]), "' target='_blank'>", result_paraloc$Gene, "</a>")),
                                   "-")
        
      }

      
      return(list("result" = result, "result_paraloc" = result_paraloc))
    } else {
      return(list("result" = result, "result_paraloc" = result_paraloc))
    }
  }
  
  # observe({
  #   if (req(input$tabs) == "Known pathogenic variants in paralogous positions")
  #     updateQueryString(paste0("?",input$mutation[1],"paralogs"), mode = "push")
  #   if (req(input$tabs) == "Paralogous Positions")
  #     updateQueryString(paste0("?",input$mutation[1],"paraloc"), mode = "push")
  # })
  
  observeEvent(input$submit_button, {
    # updateQueryString(paste0("?",input$mutation[1]), mode = "push")
    if (nrow(get_paralog("NO")$result)>=1){ # check if result table is empty
      output$paralog<-renderDataTable(DT::datatable(isolate(get_paralog("NO")$result),
                                            escape = F, # escape text hyperlink to url instead of text
                                            options = list(paging = TRUE,scrollX = TRUE, columnDefs = list(list(className = 'dt-center',targets="_all"), list(colour = 'black'))),# set options for table eg. per page lines #,columnDefs = list(list(className = 'dt-right', targets = c(1,5,7,11)))
                                            rownames = FALSE,
                                            class = "display nowrap compact",
                                            container = sketch
                                            ) %>%
                                formatStyle(c("var.query", "ID.query", "Gene.query", "Codons.query", "Transcript.query", "Protein_dot.query", "Para_Z_score.query"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold') %>%
                                formatStyle(c("Para_Z_score.query"), "border-right" = "solid 2px") %>% 
                                formatStyle(columns = colnames(.$x$data), `font-size` = "13px")
                                )
      output$paraloc<-renderDataTable(DT::datatable(isolate(get_paralog("NO")$result_paraloc),
                                                    escape = F, # escape text hyperlink to url instead of text
                                                    options = list(
                                                      searchHighlight = TRUE,
                                                      paging = TRUE,
                                                      scrollX = FALSE,
                                                      #autoWidth = TRUE,
                                                      columnDefs = list(list(width = "100px",targets = c(0)))),# set options for table eg. per page lines #,columnDefs = list(list(className = 'dt-right', targets = c(1,5,7,11)))
                                                    rownames = FALSE,
                                                    class = "display compact",
                                                    container = sketch2
                                                    ) %>%
                                        formatStyle("var", "white-space"="nowrap")
                                      )
    } else {
      #Error catching for if query returns empty table
      output$paralog<-showModal(modalDialog(
        title = "Paralog Annotator", # We can change the msg
        HTML("Your query returned no variants<br>Please try another input variant(s)<br>"), # and this msg
        easyClose = TRUE))
      shinyjs::reset("myapp") # we can delete this so the app does not restart every time
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
  output$download <- downloadHandler(
    filename = function() {
      paste("paralog_annotation", ".tsv",sep="") # need to give specific name?
    },
    content = function(file) {
      write.table(edit_output_columns(get_paralog("YES")$result), file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  
})






