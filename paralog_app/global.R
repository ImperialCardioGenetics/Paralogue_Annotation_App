library(shiny)
library(DT)
library(shinythemes)
library(stringr)

#library(tidyverse)

#read gene symbol/ENSG and write to dict
# mart_export <- read.delim(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/mart_export.txt"), quote="", stringsAsFactors=F)
mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=F)
map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)


#PRELOAD DATA ON SERVER STARTUP - THIS TAKES A WHILE - FOR TESTING BEST USE SMALLER DATASET
# raw_data = NULL
# for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
# # for (i in c(21)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
#   #use dirname(rstudioapi::getActiveDocumentContext()$path) to get relative path of this (global.R) file
#   # load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
#   # load(paste0("data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
# 
#   Total_annotations = readRDS(paste0("data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RDS"))
# 
#   if (is.null(raw_data)){
#     # Total_annotations$CHROM.x = as.character(Total_annotations$CHROM.x)
#     # Total_annotations$CHROM.y = as.character(Total_annotations$CHROM.y)
#     raw_data = Total_annotations
#   } else {
#     # Total_annotations$CHROM.x = as.character(Total_annotations$CHROM.x)
#     # Total_annotations$CHROM.y = as.character(Total_annotations$CHROM.y)
#     raw_data = base::rbind(raw_data, dplyr::setdiff(Total_annotations, raw_data))
#   }
# }
#
# saveRDS(raw_data, file = "./data/raw_data_paralog.Rds")
raw_data = readRDS("./data/raw_data_paralog.Rds")

#rm(Total_annotations)




# Paraloc_data = NULL
# for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
# # for (i in c(21)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
# 
#   Paraloc = readRDS(paste0("data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RDS"))
# 
# 
#   Paraloc = dplyr::distinct(Paraloc)
#   # Paraloc$CHROM = as.character(Paraloc$CHROM)
#   if (is.null(Paraloc_data)){
#     Paraloc_data = Paraloc
#   } else {
#     Paraloc_data = base::rbind(Paraloc_data, dplyr::setdiff(Paraloc, Paraloc_data))
#   }
# }
# 
# saveRDS(Paraloc_data, file = "./data/Paraloc_data_paraloc.Rds")
Paraloc_data = readRDS("./data/Paraloc_data_paraloc.Rds")



#rm(Paraloc)


# Paraloc_data$var = paste(Paraloc_data$CHROM,Paraloc_data$POS,Paraloc_data$REF,sep=" ")
# Paraloc_data = subset(Paraloc_data,select=c(var, Gene, Paralogue_Vars))
# Paraloc_data$Paralogue_Vars = sapply(Paraloc_data$Paralogue_Vars, stringr::str_replace, "&", "") #PROBABLY A GOOD IDEA TO DO THIS IN POST-PROCESSING BEFORE LOADING DATA IN 
# Paraloc_data$Paralogue_Vars = sapply(Paraloc_data$Paralogue_Vars, stringr::str_replace_all, "&", " ")

###LOAD DATABASE HERE
# con <- RSQLite::dbConnect(RSQLite::SQLite(), "data/db.sqlite")
###OR LOAD WHOLE RDS OBJECTS HERE
# raw_data = readRDS("data/Total_annotations_all_chrom_noQC.RDS")
# Paraloc_data = readRDS("data/Para_locations_all_chrom_noQC.RDS")


predict_output = function(input_data){
  # print(input_data)
  # print(paste(input_data$mutation, collapse = '", "'))
  # print(input_data$paraloc)

  #MOVED LOADING OF DATA TO ABOVE, OUTSIDE OF FUNCTION

  ### select the vars ## this can be done first to reduce filtering time if final dataset is huge
  ### select through SQLDB 
  # output = RSQLite::dbGetQuery(
  #   con, paste0("SELECT * FROM raw_data WHERE var IN ('",paste(input_data$mutation, collapse = "','"),"')")
  # )
  # paraloc_output = RSQLite::dbGetQuery(
  #   con, paste0("SELECT * FROM Paraloc_data WHERE var IN ('",paste(input_data$paraloc, collapse = "','"),"')")
  # )
  
  ### select through RDS objects
  output = raw_data[raw_data$var %in%  input_data$mutation,]
  paraloc_output = Paraloc_data[Paraloc_data$var %in% input_data$paraloc,]
 

  # apply protein notation change to output table only
  # output <- format_protein_notation(output, AA_map = AA_map)
  
  # will change subset to dplyr::select to change colnames at the same time
  output = dplyr::select(output,
                         var.query=var,
                         ID.query=ID,
                         Gene.query=Gene,
                         Transcript.query=Transcript,
                         Protein_dot.query=Protein_dot.x,
                         Codons.query=Codons.x,
                         Para_Z_score.query=Para_Z_score.x, 
                         var.paralog=var2, 
                         ID.paralog=ID.y, 
                         SYMBOL.paralog=SYMBOL,
                         Protein_dot.paralog=Protein_dot.y,
                         Codons.paralog=Codons.y,
                         Para_Z_score.paralog=Para_Z_score.y)
  
  #convert numeric to character so as all output df columns left align when renderDataTable()
  output = dplyr::mutate_if(output, is.numeric, as.character)
  return(list("output" = output, "paraloc_output" = paraloc_output))
  # print(output)
}

# sketch = htmltools::withTags(table(
#   class = 'display',
#   thead(
#     tr(
#       th(colspan = 7, 'Query variant(s)', 
#          # bgcolor="#cbcbcd",
#          # color = "#000000",
#          style = "border-right: solid 2px;"),
#       th(colspan = 10, 'Equivalent variant(s)')
#     ),
#     tr(
#       # lapply(c("Chrom Pos REF ALT", "ClinVar ID", "Gene", "Codons", "Protein position", "Amino acids"), th
#       lapply(c("Chrom Pos REF ALT", "ClinVar ID", "Gene", "Codons", "Transcript", "Protein"), th
#              # bgcolor="#cbcbcd", color = "#000000"
#              ),
#       th("Para_Z\nscore", style = "border-right: solid 2px;"),
#       # lapply(c("Chrom Pos REF ALT", "ClinVar ID", "Gene", "Codons", "Protein position", "Amino acids", "Para_Z score", "Ensembl alignment"), th)
#       lapply(c("Chrom Pos REF ALT", "ClinVar ID", "Gene", "Codons", "Protein", "Para_Z\nscore", "Ensembl alignment"), th)
#     )
#   )
# ))
# 
# sketch2 = htmltools::withTags(table(
#   class = 'display',
#   thead(
#     tr(
#       lapply(c("Chrom Pos REF", "Gene", "Equivalent Paralogous Locations"), th
#              # bgcolor="#cbcbcd", color = "#000000"
#       )
#     )
#   )
# ))

# use tidyr::separate to split var.query and var.paralog columns for downloaded
edit_download_cols = function(df) {
  
  df <- tidyr::separate(df, var.query, into = c("CHR.query", "POS.query", "REF.query", "ALT.query") )
  df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
  
  return(df)
}

edit_download_cols_paraloc = function(df) {
  
  df <- tidyr::separate(df, var, into = c("CHR.query", "POS.query", "REF.query") )
  #df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
  colnames(df) <- c("CHR.query", "POS.query", "REF.query", "Gene.query", "Positions.paralog")
  return(df)
}

# function to check if uploaded variants file is valid
check_upload_file = function(inFile) {
  #eg. msvcf
  #input_1row = read.table("data/head200.vcf",nrows = 1, as.is =T)
  
  #filename <- "data/test_upload2.txt"
  #filename2 <- "data/test_upload2.txt"
  
  
  # check if upload file is gz
  if (endsWith(inFile$datapath, ".gz")) {
    input_1row = read.table(gzfile(inFile$datapath),nrows = 1, as.is =T)
  } else {
    input_1row = read.table(inFile$datapath,nrows = 1, as.is =T)
  }
  
  # check if file is vcf
  if (ncol(input_1row)>=8 & grepl("^chr|^[1-9]|^[XY]", input_1row$V1) & class(input_1row$V2)=="integer" & nchar(input_1row$V1,type = "chars")<=5) {
    if (endsWith(inFile$datapath, ".gz")) {
      input_vcf = read.table(gzfile(inFile$datapath), colClasses = c("character", "integer",rep("character", 3), rep("NULL",(ncol(input_1row)-5))), stringsAsFactors = F)
    } else {
      input_vcf = read.table(inFile$datapath, colClasses = c("character", "integer",rep("character", 3), rep("NULL",(ncol(input_1row)-5))), stringsAsFactors = F)
    }
    input_file <- data.frame(do.call(paste, c(input_vcf[,1:2],input_vcf[,4:5], sep=":")),stringsAsFactors = F)
    
    # check if file is flat txt
  } else if (ncol(input_1row)>=1 & ncol(input_1row)<=4 & grepl("^chr|^[1-9]|^[XY]", input_1row$V1)) {
    if (endsWith(inFile$datapath, ".gz")) {
      input_vcf = read.table(gzfile(inFile$datapath), colClasses = "character", stringsAsFactors = F, as.is = T)
    } else {
      input_vcf = read.table(inFile$datapath, colClasses = "character",stringsAsFactors = F,as.is = T)
      #input_vcf = read.table(filename, colClasses = "character",stringsAsFactors = F,as.is = T)
      
    }
    input_file <- data.frame(do.call(paste, c(input_vcf[colnames(input_vcf)],sep=":")),stringsAsFactors = F)
    
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

# function to add ensembl URL link
add_URLs <- function(pos) {
  
  # split position string
  line <- str_split(unlist(pos[1]) , pattern = " ", simplify = T)
  # paste ensembl gene URL
  line[1,1] <- paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(line[1])]), "' target='_blank'>", unlist(line[1]), "</a>")
  # paste back to string
  pos <- paste(line[1,],collapse = " ")
  
  return(pos)
}

# function to add ensembl URLs to paralg positions in genes
add_paraloc_URL = function(result_paraloc) {
  
  # result_paraloc$Gene<- ifelse(!is.na(result_paraloc$Gene), 
  #                             (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene)]), "' target='_blank'>", result_paraloc$Gene, "</a>")),
  #                             "-")
  
  # split all positions into a vector/list
  
  #result_paraloc <- result_paraloc[complete.cases(result_paraloc),]
  
  if (nrow(result_paraloc)!=0){
    result_paraloc$Paralogue_Vars <- str_split(result_paraloc$Paralogue_Vars , pattern = ", ", simplify = F)
    
    for (i in c(1:nrow(result_paraloc))){ 
      # get all positions from list and apply add_URLs function to every position
      # then paste back as string
      result_paraloc$Paralogue_Vars[i] <- paste(unlist(lapply(unlist(result_paraloc$Paralogue_Vars[i]), function(line) add_URLs(line))), collapse = ", ")
      
    }
  } else { 
    #shiny::showNotification("No data", type = "error")
    result_paraloc <- NULL
    return(result_paraloc)    
    }

  return(result_paraloc)
}
