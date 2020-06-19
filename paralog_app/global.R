library(shiny)
library(DT)
library(shinythemes)
library(stringr)
library(tidyr)
# library(readr)
# library(tidyverse)
# library(microbenchmark)

#library(tidyverse)

#read gene symbol/ENSG and write to dict
# mart_export <- read.delim(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/mart_export.txt"), quote="", stringsAsFactors=F)
mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=F)
map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)

# Load all data as Rds data
#raw_data = readRDS("./data/raw_data_paralog.Rds")
#Paraloc_data = readRDS("../../Paraloc_data_paraloc.Rds")

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
# raw_data = readRDS("./data/raw_data_paralog.Rds")

# prepare/write files to use with tabix
# raw_data <- tidyr::separate(raw_data, var, into = c("CHR.query", "POS.query", "REF.query", "ALT.query") )
# write_tsv(raw_data, "./data/raw_data.txt")

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
# Paraloc_data = readRDS("../../Paraloc_data_paraloc.Rds")

# prepare/write files to use with tabix
# Paraloc_data <- tidyr::separate(Paraloc_data, var, into = c("CHR.query", "POS.query", "REF.query") )
# write_tsv(Paraloc_data, "./data/paraloc_data.txt")


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

#MOVED LOADING OF DATA TO ABOVE, OUTSIDE OF FUNCTION

### select the vars ## this can be done first to reduce filtering time if final dataset is huge
### select through SQLDB 
# output = RSQLite::dbGetQuery(
#   con, paste0("SELECT * FROM raw_data WHERE var IN ('",paste(input_data$mutation, collapse = "','"),"')")
# )
# paraloc_output = RSQLite::dbGetQuery(
#   con, paste0("SELECT * FROM Paraloc_data WHERE var IN ('",paste(input_data$paraloc, collapse = "','"),"')")
# )

paralog_colnames <- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  "ALT.query",
  "ID.query",
  "Gene.query",
  "Codons.query",
  "Transcript.query",
  "Protein_dot.query",
  "Para_Z_score.query",
  "var.paralog",
  "ID.paralog",
  "SYMBOL.paralog",
  "Codons.paralog",
  "Protein_dot.paralog",
  "Para_Z_score.paralog"
)
paraloc_colnames<- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  "Gene.query",
  "Positions.paralog")

lookup_paralog <- function(input_data){
  
  paralog_out <- NULL
  #paralog_out <- read.csv(text = paste(paralog_colnames,collapse = ","))
  
  
  for (i in 1:nrow(input_data)) {
    
    query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
    #CMD_paralog<- paste0("tabix ", paralog_data, " ", query)
    #tabix_paralog <- system(command = paste0("tabix ", paralog_data, " ", query), intern = T,wait = T)
    
    # maybe need to set tbaxi dir
    # tabix <- paste("")
    tabix_paralog <- system(command = paste0("tabix data/raw_data_sorted.txt.gz ", query), intern = T,wait = T)
    
    #pg1 <- separate(as.data.frame(unlist(tabix_paralog)),1,sep = "\t", into =  paralog_colnames)
    #pg1 <- pg1[(pg1$REF.query==input_data$REF.query[i] & pg1$ALT.query == input_data$ALT.query[i]),]
    
    
    pg1 <- as.data.frame(stringr::str_split_fixed(tabix_paralog, pattern = "\t", length(paralog_colnames)), stringsAsFactors = F)
    #colnames(pg1) <- paralog_colnames
    
    
    if (is.null(paralog_out)){
      paralog_out = pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),]
    } else {
      paralog_out = rbind(paralog_out, pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),])
    }
  }
  #paralog_out = rbind(paralog_out, pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),])
  
  colnames(paralog_out) <- paralog_colnames
  #paralog_out <- dplyr::na_if(paralog_out, "NA")
  
  
  # if (nrow(paralog_out)!=0) {
  #   paralog_out <- dplyr::na_if(paralog_out, "NA")
  #   return(paralog_out)
  # } else {
  #   paralog_out <- NULL
  #   return(paralog_out)
  # }
  return(paralog_out)
}

# call lookup function 
# paralog_out<- lookup_paralog(input_data)


lookup_paraloc <- function(input_data){
  
  paraloc_out <- NULL
  
  for (i in 1:nrow(input_data)) {
    
    query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
    #CMD_paraloc<- paste0("tabix ", paraloc_data, " ", query)
    #tabix_paraloc <- system(command =  paste0("tabix ", paraloc_data, " ", query), intern = T, wait = T)
    tabix_paraloc <- system(command =  paste0("tabix data/paraloc_data_sorted_fix_header.txt.gz " , query), intern = T, wait = T)
    
    #pc1 <- as.data.frame(t(unlist(strsplit(tabix_paraloc,split="\t"))),stringsAsFactors = F)
    pc1<- as.data.frame(stringr::str_split_fixed(tabix_paraloc, pattern = "\t", n = length(paraloc_colnames)), stringsAsFactors = F)
    #colnames(pc1) <- paraloc_colnames
    
    #pc1 <- pc1[(pc1$REF.query==input_data$REF.query[i]),]
    
    if (is.null(paraloc_out)){
      paraloc_out = pc1[(pc1$V3==input_data$REF.query[i]),]
    } else {
      paraloc_out = rbind(paraloc_out, pc1[(pc1$V3==input_data$REF.query[i]),])
    }
  }
  #paraloc_out = rbind(paraloc_out, pc1[(pc1$V3==input_data$REF.query[i]),])
  
  colnames(paraloc_out) <- paraloc_colnames
  # if (nrow(paraloc_out)!=0) {
  #   paraloc_out <- dplyr::na_if(paraloc_out, "NA")
  #   return(paraloc_out)
  # } else {
  #   return(paraloc_out)
  # }
  return(paraloc_out)
  
}




predict_out = function(input_data){


  #paralog_out <- read.csv(text = paste(paralog_colnames,collapse = ","))

  input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F) 

  
  
  # call lookup function 
  # paraloc_out <- lookup_paraloc(input_data)
  
  return(list("output" = lookup_paralog(input_data), "paraloc_output" = lookup_paraloc(input_data)))
  
}  


#result <- predict_out(input_data = input_data)

  
predict_output = function(input_data){
  #### select through RDS objects ####
  paralog_output = raw_data[raw_data$var %in%  input_data$mutation,]
  paraloc_output = Paraloc_data[Paraloc_data$var %in% input_data$paraloc,]
  
  # will change subset to dplyr::select to change colnames at the same time
  paralog_output = dplyr::select(paralog_output,
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
  paralog_output = dplyr::mutate_if(paralog_output, is.numeric, as.character)
  
  return(list("output" = paralog_output, "paraloc_output" = paraloc_output))
  # print(output)
}

#result <- predict_output(input_data = input_data)



# # BENCHMARK
# mb2 <- microbenchmark(out_new = predict_out(input_data),
#                      output_old = predict_output(input_data), times = 10)
# 
# 
# mb_pg <- microbenchmark(
#   pg1 = separate(as.data.frame(unlist(tabix_paralog)),1,sep = "\t", into =  paralog_colnames),
#   pg2 = as.data.frame(str_split_fixed(tabix_paralog, pattern = "\t", n = length(paralog_colnames)), stringsAsFactors = F),
# times = 30)
# 
# mb <- microbenchmark(paralog_out = lookup_paralog(input_data),
#                       paraloc_out = lookup_paraloc(input_data),
#                       paralog_output = raw_data[raw_data$var %in%  input_data$mutation,],
#                       paraloc_output = Paraloc_data[Paraloc_data$var %in% input_data$paraloc,], times = 10)
# 
# 
# autoplot(mb)
# autoplot(mb_pg)
# autoplot(mb2)
# 


# use tidyr::separate to split var.query and var.paralog columns for downloaded
edit_download_cols = function(df) {
  
  #df <- tidyr::separate(df, var.query, into = c("CHR.query", "POS.query", "REF.query", "ALT.query") )
  df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
  
  return(df)
}

# edit_download_cols_paraloc = function(df) {
#   
#   df <- tidyr::separate(df, var, into = c("CHR.query", "POS.query", "REF.query") )
#   #df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
#   colnames(df) <- c("CHR.query", "POS.query", "REF.query", "Gene.query", "Positions.paralog")
#   return(df)
# }

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
    result_paraloc$Positions.paralog <- str_split(result_paraloc$Positions.paralog , pattern = ", ", simplify = F)
    
    for (i in c(1:nrow(result_paraloc))){ 
      # get all positions from list and apply add_URLs function to every position
      # then paste back as string
      result_paraloc$Positions.paralog[i] <- paste(unlist(lapply(unlist(result_paraloc$Positions.paralog[i]), function(line) add_URLs(line))), collapse = ", ")
      
    }
  } else { 
    #shiny::showNotification("No data", type = "error")
    result_paraloc <- NULL
    return(result_paraloc)    
    }

  return(result_paraloc)
}
