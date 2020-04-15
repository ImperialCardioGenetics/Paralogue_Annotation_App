library(shiny)
library(DT)
library(shinythemes)

#library(tidyverse)

#PRELOAD DATA ON SERVER STARTUP - THIS TAKES A WHILE - FOR TESTING BEST USE SMALLER DATASET
raw_data = NULL
# for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
for (i in c(1)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
  load(paste0("data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
  if (is.null(raw_data)){
    Total_annotations$CHROM.x = as.character(Total_annotations$CHROM.x)
    Total_annotations$CHROM.y = as.character(Total_annotations$CHROM.y)
    raw_data = Total_annotations
  } else {
    Total_annotations$CHROM.x = as.character(Total_annotations$CHROM.x)
    Total_annotations$CHROM.y = as.character(Total_annotations$CHROM.y)
    raw_data = base::rbind(raw_data, dplyr::setdiff(Total_annotations, raw_data))
  }
}
raw_data$var = paste(raw_data$CHROM.x,raw_data$POS.x,raw_data$REF.x,raw_data$ALT.x,sep="\t")
raw_data = subset(raw_data,select=c(var, Gene, Codons.x, Protein_position.x, Amino_acids.x, Para_Z_score.x, CHROM.y, POS.y, REF.y, ALT.y, ID.y, SYMBOL, Codons.y, Protein_position.y, Amino_acids.y, Para_Z_score.y))


predict_output = function(input_data){
  print(input_data)
  # I opened formated and outputted again the RData odject from Nick to edit rownames and NAs 
  # This can be done here within a function
  # load the RData or rds file
  # load("data/paralog_tmp.RData")
  
  #   exported RData with hardcoded dataframe name raw_data
  # load("data/place_holder_results_datatable.RData",)
  
  #MOVED LOADING OF DATA TO ABOVE, OUTSIDE OF FUNCTION

  # write new column chr:pos:ref:alt to look up
  # raw_data$var<-paste(raw_data$CHROM.x,raw_data$POS.x,raw_data$REF.x,raw_data$ALT.x,sep="\t") #done above as website loads
  
  #generate Clinvar URLs for query and result
  # moved that to the server side
  #raw_data$link.x<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",raw_data$ID.x,"/"), "' target='_blank'>", raw_data$ID.x, "</a>")  
  #raw_data$link.x<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",raw_data$ID.y,"/"), "' target='_blank'>", raw_data$ID.y, "</a>")  
  
  # select dataframe columns if not formated
  # paralog_tmp<-subset(raw_data,select=c(raw_data$var,raw_data$CHROM.y ,raw_data$POS.y,raw_data$REF.y,raw_data$ALT.y,raw_data$ID.y ,raw_data$SYMBOL,raw_data$Protein_position.y,raw_data$REF_Amino_acids.y,raw_data$ALT_Amino_acids.y ,raw_data$Codons.y,raw_data$Para_Z_score.y))
  # raw_data<-subset(raw_data,select=c(var, Gene, ID.x, CHROM.y ,POS.y,REF.y,ALT.y,ID.y ,SYMBOL,Protein_position.y,REF_Amino_acids.y,ALT_Amino_acids.y ,Codons.y,Para_Z_score.y),)
  
  # rename dataframe columns for webpage
  # colnames(raw_data)<-c("Variant ID","Query_Gene","Query_ClinVar", "Chr","Position","REF","ALT","ClinVar_ID","Gene","Protein Position","Reference AA", "Alt AA","Codons","para_Z Score" )

  # select the vars ## this can be done first to reduce filtering time if final dataset is huge
  output = raw_data[raw_data$var %in% input_data$mutation,]
  # output = raw_data[raw_data$CHROM.x == as.numeric(input_data$chr) & 
  #                     raw_data$POS.x == as.numeric(input_data$pos) & 
  #                     raw_data$REF.x == input_data$ref &
  #                     raw_data$ALT.x == input_data$alt,]
  output$CHROM.x = sapply(strsplit(output$var, "\t"), "[", 1)
  output$POS.x = sapply(strsplit(output$var, "\t"), "[", 2)
  output$REF.x = sapply(strsplit(output$var, "\t"), "[", 3)
  output$ALT.x = sapply(strsplit(output$var, "\t"), "[", 4)
  output = subset(output,select=c(CHROM.x, POS.x, REF.x, ALT.x, Gene, Codons.x, Protein_position.x, Amino_acids.x, Para_Z_score.x, CHROM.y, POS.y, REF.y, ALT.y, ID.y, SYMBOL, Codons.y, Protein_position.y, Amino_acids.y, Para_Z_score.y))
  
  #in order to remove duplicated query rows
  # tmp_chrom = NULL
  # tmp_pos = NULL
  # tmp_ref = NULL
  # tmp_alt = NULL
  # for (row in 1:nrow(output)){
  #   
  # }
  
  return(output)
  print(output)
}

#use dirname(rstudioapi::getActiveDocumentContext()$path) to get relative path of this (global.R) file
clinvar_P_LP = read.csv(file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/clinvar/clinvar_20190114_GRCh37_onlyPathogenic_and_Likely_pathogenic.vcf"), sep = "\t", comment.char = "#", stringsAsFactors = F, header = F) #load in clinvar data for query variant
clinvar_P_LP = clinvar_P_LP[,1:5]
colnames(clinvar_P_LP) = c("CHR", "POS", "ID", "REF", "ALT")
clinvar_P_LP$var = paste(clinvar_P_LP$CHR,clinvar_P_LP$POS,clinvar_P_LP$REF,clinvar_P_LP$ALT,sep="\t")
clinvar_P_LP = subset(clinvar_P_LP,select=c(var, ID))

predict_output_for_known = function(input_data){
  #print(input_data)
  output = clinvar_P_LP[clinvar_P_LP$var %in% input_data$mutation,]
  output$CHR = sapply(strsplit(output$var, "\t"), "[", 1)
  output$POS = sapply(strsplit(output$var, "\t"), "[", 2)
  output$REF = sapply(strsplit(output$var, "\t"), "[", 3)
  output$ALT = sapply(strsplit(output$var, "\t"), "[", 4)
  output = subset(output,select=c(CHR, POS, ID, REF, ALT))
  
  output$ID<- paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",output$ID,"/"), "' target='_blank'>", output$ID, "</a>")  
  
  print(output)
  return(output)
}

#MAY NOT EVEN NEED THIS FUNCTION BELOW, COULD JUST INTEGRATE TO FUNCTION ABOVE
check_if_known = function(chr,pos,ref,alt){
  query_variant_ID = clinvar_P_LP[clinvar_P_LP$CHR == chr & clinvar_P_LP$POS == pos & clinvar_P_LP$REF == ref & clinvar_P_LP$ALT == alt,]
  if (length(query_variant_ID)>0){
    query_variant_ID = paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",query_variant_ID,"/"), "' target='_blank'>", query_variant_ID, "</a>")  
  } else {
    query_variant_ID = "Th"
  }
  return(query_variant_ID)
}



sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(colspan = 9, 'Query variant(s)', 
         # bgcolor="#cbcbcd",
         # color = "#000000",
         style = "border-right: solid 2px;"),
      th(colspan = 11, 'Equivalent variant(s)')
    ),
    tr(
      lapply(c("Chr", "Position", "REF", "ALT", "Gene", "Codons", "Protein position", "Amino acids"), th
             # bgcolor="#cbcbcd", color = "#000000"
             ),
      th("Para_Z score", style = "border-right: solid 2px;"),
      lapply(c("Chr", "Position", "REF", "ALT", "ClinVar ID", "Gene", "Codons", "Protein position", "Amino acids", "Para_Z score", "Ensembl alignment"), th)
    )
  )
))

sketch2 = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      lapply(c("Chr", "Position", "ClinVar ID", "REF", "ALT"), th
             # bgcolor="#cbcbcd", color = "#000000"
      )
    )
  )
))
