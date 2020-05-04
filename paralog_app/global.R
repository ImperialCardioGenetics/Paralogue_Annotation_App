library(shiny)
library(DT)
library(shinythemes)
library(stringr)

#library(tidyverse)

#read gene symbol/ENSG and write to dict
mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=FALSE)
map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)


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
  #Paraloc$var = paste(Paraloc$CHROM,Paraloc$POS,Paraloc$REF,Paraloc$Gene,sep=" ")
  #Paraloc = subset(Paraloc,select=c(var, Paralogue_Vars))
  Paraloc = subset(Paraloc, select=c(CHROM,POS,REF,Gene,Paralogue_Vars)) #IF NOT COMBINING INTO VAR THEN NEED TO CHANGE HOW WE LOOK UP DATA
  Paraloc = dplyr::distinct(Paraloc)
  if (is.null(Paraloc_data)){
    Paraloc_data = Paraloc
  } else {
    Paraloc_data = base::rbind(Paraloc_data, dplyr::setdiff(Paraloc, Paraloc_data))
  }
}
rm(Paraloc)
# Paraloc_data$var = paste(Paraloc$CHROM,Paraloc$POS,Paraloc$REF,Paraloc$ALT,sep=" ")
# Paraloc_data = subset(Paraloc_data,select=c(var, Gene, Paralogue_Vars))
# Paraloc_data$Paralogue_Vars = sapply(Paraloc_data$Paralogue_Vars, stringr::str_replace, "&", "") #PROBABLY A GOOD IDEA TO DO THIS IN POST-PROCESSING BEFORE LOADING DATA IN 
# Paraloc_data$Paralogue_Vars = sapply(Paraloc_data$Paralogue_Vars, stringr::str_replace_all, "&", " ")

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

# function to check if uploaded variants file is valid
check_upload_file = function(inFile) {
  #eg. msvcf
  #input_1row = read.table("data/head200.vcf",nrows = 1, as.is =T)
  
  # check if upload file is gz
  if (endsWith(inFile$datapath, ".gz")) {
    input_1row = read.table(gzfile(inFile$datapath),nrows = 1, as.is =T)
  } else {
    input_1row = read.table(inFile$datapath,nrows = 1, as.is =T)
  }
  
  # check if file is vcf
  if (ncol(input_1row)>=8 & grepl("^chr|^[1-9]|^[XY]", input_1row$V1) & class(input_1row$V2)=="integer" & nchar(input_1row$V1,type = "chars")<=5) {
    if (endsWith(inFile$datapath, ".gz")) {
      input_vcf = read.table(gzfile(inFile$datapath), colClasses = c("character", "integer",rep("character", 3), rep("NULL",(ncol(input_1row)-5))))
    } else {
      input_vcf = read.table(inFile$datapath, colClasses = c("character", "integer",rep("character", 3), rep("NULL",(ncol(input_1row)-5))))
    }
    input_file <- data.frame(do.call(paste, c(input_vcf[,1:2],input_vcf[,4:5], sep=" ")),stringsAsFactors = F)
    
    # check if file is flat txt
  } else if (ncol(input_1row)>=1 & ncol(input_1row)<=4 & grepl("^chr|^[1-9]|^[XY]", input_1row$V1)) {
    if (endsWith(inFile$datapath, ".gz")) {
      input_vcf = read.table(gzfile(inFile$datapath), colClasses = "character")
    } else {
      input_vcf = read.table(inFile$datapath, colClasses = "character")
    }
    input_file <- data.frame(do.call(paste, c(input_vcf[colnames(input_vcf)],sep=" ")),stringsAsFactors = F)
    
    # check if table if correct format
    if (nchar(input_file[1,1],type = "chars")<10 | !grepl("[AGTC]$", input_file[1,1])) {
      
      # write NA table
      input_file<- data.frame(col1=NA)
    }
  } else {
    # write NA table
    input_file<- data.frame(col1=NA)
  }
  
  return(input_file)
}

