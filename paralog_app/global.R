library(shiny)
library(DT)
library(shinythemes)

#library(tidyverse)

#PRELOAD DATA ON SERVER STARTUP - THIS TAKES A WHILE - FOR TESTING BEST USE SMALLER DATASET
raw_data = NULL
# for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
for (i in c(1)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
  #use dirname(rstudioapi::getActiveDocumentContext()$path) to get relative path of this (global.R) file
  load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
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
raw_data$var = paste(raw_data$CHROM.x,raw_data$POS.x,raw_data$REF.x,raw_data$ALT.x,sep=" ")
raw_data$var2 = paste(raw_data$CHROM.y,raw_data$POS.y,raw_data$REF.y,raw_data$ALT.y,sep=" ")
raw_data = subset(raw_data,select=c(var, Gene, Codons.x, Protein_position.x, Amino_acids.x, Para_Z_score.x, var2, ID.y, SYMBOL, Codons.y, Protein_position.y, Amino_acids.y, Para_Z_score.y))

#use dirname(rstudioapi::getActiveDocumentContext()$path) to get relative path of this (global.R) file
clinvar_P_LP = read.csv(file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/clinvar/clinvar_20190114_GRCh37_onlyPathogenic_and_Likely_pathogenic.vcf"), sep = "\t", comment.char = "#", stringsAsFactors = F, header = F) #load in clinvar data for query variant
clinvar_P_LP = clinvar_P_LP[,1:5]
colnames(clinvar_P_LP) = c("CHR", "POS", "ID", "REF", "ALT")
clinvar_P_LP$var = paste(clinvar_P_LP$CHR,clinvar_P_LP$POS,clinvar_P_LP$REF,clinvar_P_LP$ALT,sep=" ")
clinvar_P_LP = subset(clinvar_P_LP,select=c(var, ID))

raw_data = dplyr::right_join(clinvar_P_LP,raw_data,by = c("var"))

Paraloc_data = NULL
# for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
for (i in c(1)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
  #use dirname(rstudioapi::getActiveDocumentContext()$path) to get relative path of this (global.R) file
  load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
  if (is.null(Paraloc_data)){
    Paraloc_data = Paraloc
  } else {
    Paraloc_data = base::rbind(Paraloc_data, dplyr::setdiff(Paraloc, Paraloc_data))
  }
}
#rm(Paraloc)
Paraloc_data$var = paste(Paraloc$CHROM,Paraloc$POS,Paraloc$REF,Paraloc$ALT,sep=" ")
Paraloc_data = subset(Paraloc_data,select=c(var, Gene, Paralogue_Vars))

predict_output = function(input_data){
  print(input_data$mutation)
  print(raw_data$var[1])
  
  #MOVED LOADING OF DATA TO ABOVE, OUTSIDE OF FUNCTION

  # select the vars ## this can be done first to reduce filtering time if final dataset is huge
  output = raw_data[raw_data$var %in%  input_data$mutation,]
  paraloc_output = Paraloc_data[Paraloc_data$var %in% input_data$mutation,]
  
  # will change subset to dplyr::select to change colnames at the same time
  output = dplyr::select(output,
    var.query=var,
    ID.query=ID,
    Gene.query=Gene, 
    Codons.query=Codons.x, 
    Protein_position.query=Protein_position.x, 
    Amino_acids.query=Amino_acids.x, 
    Para_Z_score.query=Para_Z_score.x, 
    var.paralog=var2, 
    ID.paralog=ID.y, 
    SYMBOL.paralog=SYMBOL, 
    Codons.paralog=Codons.y, 
    Protein_position.paralog=Protein_position.y, 
    Amino_acids.paralog=Amino_acids.y, 
    Para_Z_score.paralog=Para_Z_score.y)
  
  #convert numeric to character so as all output df columns left align when renderDataTable()
  output = dplyr::mutate_if(output, is.numeric, as.character)
  return(list("output" = output, "paraloc_output" = paraloc_output))
  print(output)
}

sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(colspan = 7, 'Query variant(s)', 
         # bgcolor="#cbcbcd",
         # color = "#000000",
         style = "border-right: solid 2px;"),
      th(colspan = 11, 'Equivalent variant(s)')
    ),
    tr(
      lapply(c("Chrom Pos REF ALT", "ClinVar ID", "Gene", "Codons", "Protein position", "Amino acids"), th
             # bgcolor="#cbcbcd", color = "#000000"
             ),
      th("Para_Z score", style = "border-right: solid 2px;"),
      lapply(c("Chrom Pos REF ALT", "ClinVar ID", "Gene", "Codons", "Protein position", "Amino acids", "Para_Z score", "Ensembl alignment"), th)
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

# use tidyr::separate to split var.query and var.paralog columns for downloaded
edit_output_columns = function(df) {
  
  df <- tidyr::separate(df, var.query, into = c("CHR.query", "POS.query", "REF.query", "ALT.query") )
  df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
  
  return(df)
}





