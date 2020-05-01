library(shiny)
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)

#library(tidyverse)
options(shiny.maxRequestSize=200*1024^2) #max upload size = 200 mb

shinyServer(function(input, output){
  
  #read gene symbol/ENSG and write to dict
  mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=FALSE)
  map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)
  
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
        input_data$mutation<-gsub(":"," ",input_data$mutation)
        input_data$mutation<-gsub("^chr","",input_data$mutation)
        colnames(input_data)<-"mutation"
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
        result$ID.paralog<- ifelse(!is.na(result$ID.paralog), (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.paralog,"/"), "' target='_blank'>", result$ID.paralog, "</a>")) , NA)
        
        #ClinVarID query URL
        result$ID.query<- ifelse(!is.na(result$ID.query), (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.query,"/"), "' target='_blank'>", result$ID.query, "</a>")) , NA)
        
        #Ensembl alignment URL
        # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
        result$Ensembl_alignment_link<- ifelse(!is.na(result$SYMBOL), (paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene.query)],";g1=",map[unlist(result$SYMBOL.paralog)]), "' target='_blank'>alignment</a>")) , NA) 
      }

      
      return(list("result" = result, "result_paraloc" = result_paraloc))
    } else {
      return(list("result" = result, "result_paraloc" = result_paraloc))
    }
  }
  
  observeEvent(input$submit_button, {
    if (nrow(get_paralog("NO")$result)>=1){ # check if result table is empty
      
      output$paralog<-renderDataTable(DT::datatable(isolate(get_paralog("NO")$result),
                                            escape = F, # escape text hyperlink to url instead of text
                                            options = list(paging = TRUE,scrollX = TRUE),# set options for table eg. per page lines
                                            rownames = FALSE,
                                            class = "display nowrap",
                                            container = sketch
                                            ) %>%
                                formatStyle(c("var.query", "ID.query", "Gene.query", "Codons.query", "Protein_position.query", "Amino_acids.query", "Para_Z_score.query"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold') %>%
                                formatStyle(c("Para_Z_score.query"), "border-right" = "solid 2px")
                                )
      output$paraloc<-renderDataTable(DT::datatable(isolate(get_paralog("NO")$result_paraloc)
                                                    )
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






