library(shiny)
library(DT)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)

#library(tidyverse)

shinyServer(function(input, output){
  
  #read gene symbol/ENSG and write to dict
  mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=FALSE)
  map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)
  
  get_paralog<-function(savefile){
    
    #input<-data.frame(chr="1",pos="114713907",ref="T",alt="A")
    if(input$format=='pick'){
      req(input$chr)
      req(input$pos)
      req(input$ref)
      req(input$alt)
      var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
      var = var[nzchar(x=var)]
      input_data = data.frame(mutation=var, stringsAsFactors = FALSE)
      # input_data = data.frame(chr = input$chr, pos = input$pos, ref = input$ref, alt = input$alt) #not needed
      result<-predict_output(input_data)
    }else{
      if(input$format=='paste'){
      #input<-data.frame(var="1:114713907:T:G",stringsAsFactors = F)  
      #input$var<-data.frame(var="1\t114713907\tT\tG")
        req(input$var)
        # print(input$var)
        var<-unlist(strsplit(input$var,split="\n"))
        # print(var)
        var=var[nzchar(x=var)]
        input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
        input_data$mutation<-gsub(":"," ",input_data$mutation)
        colnames(input_data)<-"mutation"
        result<-predict_output(input_data)
    }else{
      if(input$format == 'upload') {
        req(input$file)
        inFile <- input$file
        input_file = read.table(inFile$datapath)
        input_file$V1<-gsub(":"," ",input_file$V1)
        colnames(input_file) <- "mutation"
        result <- predict_output(input_file)
      }
    }
  }
    if (savefile=="NO"){
      #add ClinVar IDs with URLs 
      #result$Query_ClinVar_link<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$Query_ClinVar,"/"), "' target='_blank'>", result$Query_ClinVar, "</a>")  #Not possible for custom_ids; Can add feature to check P/LP tableized file
      if (nrow(result)!=0){ # that where the error was generated
        result$ID.y<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.y,"/"), "' target='_blank'>", result$ID.y, "</a>")  
      
      #generate ensembl alignment URLs
      # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
        result$Ensembl_alignment_link<- paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene)],";g1=",map[unlist(result$SYMBOL)]), "' target='_blank'>alignment</a>")  
      }
      #edit and add this later
      # reorder_cols<-c("Variant_ID","Query_Gene","Query_ClinVar_link", "Chr","Position","REF","ALT","ClinVar_ID_link","Gene","Protein Position","Reference AA", "Alt AA","Codons","para_Z Score","Ensembl_alignment_link" )
      # rename_cols<-c("Query Variant","Query Gene","Query ClinVar ID", "Chr","Position","REF","ALT","ClinVar ID","Gene","Protein Position","REF AA", "ALT AA","Codons","para_Z Score","Ensembl alignment" )
      # result<- result[,reorder_cols]
      # names(result)<-rename_cols
      
      return(result)
    } else {
      return(result)
    }
  }
  
  observeEvent(input$sumbit_button, {
    if (nrow(get_paralog("NO"))>=1){ # check if result table is empty
      
      output$paralog<-renderDataTable(DT::datatable(isolate(get_paralog("NO")),
                                            escape = F, # escape text hyperlink to url instead of text
                                            options = list(paging = TRUE,scrollX = TRUE),# set options for table eg. per page lines
                                            rownames = FALSE,
                                            class = "display nowrap",
                                            container = sketch
                                            ) %>%
                                formatStyle(c("var", "ID", "Gene", "Codons.x", "Protein_position.x", "Amino_acids.x", "Para_Z_score.x"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold') %>%
                                formatStyle(c("Para_Z_score.x"), "border-right" = "solid 2px")
                                )
      
      # output$known_clinvar<-renderDataTable(DT::datatable(isolate(get_known()),
      #                                                     escape = F,
      #                                                     options = list(paging = FALSE,scrollX = TRUE),
      #                                                     rownames = FALSE, 
      #                                                     container = sketch2
      # ) %>%
      #   formatStyle(c("CHR", "POS", "ID", "REF", "ALT"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold')
      # )
      
      # output$paralog<-renderDataTable(DT::datatable(isolate(get_paralog("NO")),
      #                                               escape = F, # escape text hyperlink to url instead of text
      #                                               options = list(paging = FALSE,scrollX = TRUE),# set options for table eg. per page lines
      #                                               rownames = FALSE,
      #                                               container = sketch,
      #                                               caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;','Table 1 : ', htmltools::em('Paralogous Variants'))
      # ) %>%
      #   formatStyle(c("CHROM.x", "POS.x", "REF.x", "ALT.x", "Gene", "Codons.x", "Protein_position.x", "Amino_acids.x", "Para_Z_score.x"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold') %>%
      #   formatStyle(c("Para_Z_score.x"), "border-right" = "solid 2px")
      # )
      
    } else {
      #ERROR CATCHING FOR IF QUERY IS KNOWN
      output$known_clinvar<- showModal(modalDialog(
        title = "The input variant(s) were not found in ClinVar", # We can change the msg
        "Please try another input variant(s)", # and this msg
        easyClose = TRUE))
      shinyjs::reset("myapp") # we can delete this so the app does not restart every time
      
      #ERROR CATCHING FOR IF QUERY HAS EQUIVALENT PARALOGUE LOCATION AND/OR VARIANT 
      output$paralog<- showModal(modalDialog(
        title = "No paralogous variants found", # We can change the msg
        "Please try another input variant(s)", # and this msg
        easyClose = TRUE))
      shinyjs::reset("myapp") # we can delete this so the app does not restart every time
      
    }
    
  })
  
  observe({
    if (is.null(input$var) || input$var == "") {
      shinyjs::disable("sumbit_button")
      shinyjs::disable("reset_button")
    } else {
      shinyjs::enable("sumbit_button")
      shinyjs::enable("reset_button")
    }
  })
  
  observeEvent(input$reset, {
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
      write.table(get_paralog("YES"), file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  
})






