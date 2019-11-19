library(shiny)
library(DT)
library(shinythemes)

#library(tidyverse)

shinyServer(function(input, output){
  
  #read gene symbol/ENSG and write to dict
  mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=FALSE)
  map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)
  
  get_paralog<-function(){
    
    #input<-data.frame(chr="1",pos="114713907",ref="T",alt="A")
    
    if(input$format=='pick'){
      req(input$chr)
      req(input$pos)
      req(input$ref)
      req(input$alt)
      var<-paste(input$chr,input$pos,input$ref,input$alt,sep = ":")
      var=var[nzchar(x=var)]
      input_data<-data.frame(mutation=var)
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
  
    #add ClinVar IDs with URLs 
    result$Query_ClinVar<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$Query_ClinVar,"/"), "' target='_blank'>", result$Query_ClinVar, "</a>")  
    result$ClinVar_ID<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ClinVar_ID,"/"), "' target='_blank'>", result$ClinVar_ID, "</a>")  
    
    #generate ensembl alignment URLs
    # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
    result$Ensembl_alignment<- paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Query_Gene)],";g1=",map[unlist(result$Gene)]), "' target='_blank'>alignment</a>")  
  
    return(result)
  }
  
  output$paralog<-renderDataTable(DT::datatable(get_paralog(),
                                                escape = F, # escape text hyperlink to url instead of text
                                                options = list(paging = FALSE),# set options for table eg. per page lines
                                                rownames = F, 
                                                caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;','Table 1 : ', htmltools::em('Paralogous Variants'))
                                                ) %>%
                                                formatStyle(c("Variant_ID","Query_Gene","Query_ClinVar"),  color = 'black', backgroundColor = 'lightgrey', fontWeight = 'bold')
                                  )
  
  output$download <- downloadHandler(
    filename = function() {
      paste("paralog_annotation", ".tsv",sep="") # need to give specific name?
    },
    content = function(file) {
      write.table(get_paralog(), file, row.names = FALSE,quote = F,sep="\t")
    }
  )
  
})






