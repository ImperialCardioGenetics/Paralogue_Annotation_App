library(shiny)
library(DT)
library(shinythemes)

#library(tidyverse)



predict_output<-function(output,input_data){
  
  # I opened formated and outputted again the RData odject from Nick to edit rownames and NAs 
  # This can be done here within a function
  # load the RData or rds file
  # load("data/paralog_tmp.RData")
  
  #   exported RData with hardcoded dataframe name raw_data
  load("data/place_holder_results_datatable.RData",)
  raw_data<-place_holder_results_datatable
  
  # write new column chr:pos:ref:alt to look up
  raw_data$var<-paste(raw_data$CHROM.x,raw_data$POS.x,raw_data$REF.x,raw_data$ALT.x,sep=":")
  
  #generate Clinvar URLs for query and result
  # moved that to the server side
  #raw_data$link.x<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",raw_data$ID.x,"/"), "' target='_blank'>", raw_data$ID.x, "</a>")  
  #raw_data$link.x<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",raw_data$ID.y,"/"), "' target='_blank'>", raw_data$ID.y, "</a>")  
  
  # select dataframe columns if not formated
  # paralog_tmp<-subset(raw_data,select=c(raw_data$var,raw_data$CHROM.y ,raw_data$POS.y,raw_data$REF.y,raw_data$ALT.y,raw_data$ID.y ,raw_data$SYMBOL,raw_data$Protein_position.y,raw_data$REF_Amino_acids.y,raw_data$ALT_Amino_acids.y ,raw_data$Codons.y,raw_data$Para_Z_score.y))
  raw_data<-subset(raw_data,select=c(var, Gene, ID.x, CHROM.y ,POS.y,REF.y,ALT.y,ID.y ,SYMBOL,Protein_position.y,
                                     REF_Amino_acids.y,ALT_Amino_acids.y ,Codons.y,Para_Z_score.y),)
  # rename dataframe columns for webpage
  colnames(raw_data)<-c("Variant_ID","Query_Gene","Query_ClinVar", "Chr","Position","REF","ALT","ClinVar_ID","Gene","Protein Position","Reference AA", "Alt AA","Codons","para_Z Score" )
  
  # select the vars ## this can be done first to reduce filtering time if final dataset is huge
  output<- raw_data[raw_data$Variant_ID %in% input_data$mutation,]
  
  return(output)
}

