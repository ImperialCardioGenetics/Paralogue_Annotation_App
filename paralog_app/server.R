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
      var = paste(input$chr,input$pos,input$ref,input$alt,sep = "\t")
      var = var[nzchar(x=var)]
      input_data = data.frame(mutation=var, stringsAsFactors = FALSE)
      # input_data = data.frame(chr = input$chr, pos = input$pos, ref = input$ref, alt = input$alt) #not needed
      result<-predict_output(raw_data,input_data)
  }else{
    if(input$format=='paste'){
      req(input$var)
      var<-unlist(strsplit(input$var,split="\\s+"))
      var=var[nzchar(x=var)]
      input_data<-data.frame(mutation=var)
      colnames(input_data)<-"mutation"
      result<-predict_output(raw_data,input_data)
    }else{
        if (input$format == 'upload') {
          req(input$file)
          inFile <- input$file
          input_file = read.table(inFile$datapath)
          colnames(input_file) <- "mutation"
          result <- predict_output(raw_data, input_file)
        }
      }
    }
    if (savefile=="NO"){
      #add ClinVar IDs with URLs 
      #result$Query_ClinVar_link<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$Query_ClinVar,"/"), "' target='_blank'>", result$Query_ClinVar, "</a>")  #Not possible for custom_ids; Can add feature to check P/LP tableized file
      result$ID.y<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.y,"/"), "' target='_blank'>", result$ID.y, "</a>")  
      
      #generate ensembl alignment URLs
      # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
      result$Ensembl_alignment_link<- paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene)],";g1=",map[unlist(result$SYMBOL)]), "' target='_blank'>alignment</a>")  
      
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
    output$paralog<-renderDataTable(DT::datatable(isolate(get_paralog("NO")),
                                                escape = F, # escape text hyperlink to url instead of text
                                                options = list(paging = FALSE,scrollX = TRUE),# set options for table eg. per page lines
                                                rownames = FALSE, 
                                                container = sketch,
                                                caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;','Table 1 : ', htmltools::em('Paralogous Variants'))
                                                ) %>%
                                                formatStyle(c("CHROM.x", "POS.x", "REF.x", "ALT.x", "Gene", "Codons.x", "Protein_position.x", "Amino_acids.x", "Para_Z_score.x"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold') %>%
                                      formatStyle(c("Para_Z_score.x"), "border-right" = "solid 2px")
                                  )
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






