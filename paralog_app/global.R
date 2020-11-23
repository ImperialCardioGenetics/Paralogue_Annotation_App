library(shiny)
library(DT)
library(shinythemes)
library(tidyverse)
#library(microbenchmark)
#library(drawProteins)
library(plotly)
library(grid)
library(httr)





#read gene symbol/ENSG and write to dict
# mart_export <- read.delim(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/mart_export.txt"), quote="", stringsAsFactors=F)
mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=F)
map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)

HGNC_export <- read.delim("data/HGNC_all_genes.txt", quote="", stringsAsFactors=F)
# HGNC_export <- tidyr::separate(HGNC_export,HGNC_ID, into = c("HGNC", "ID"), remove = T)
map_HGNC = setNames(HGNC_export$UniProt_ID, HGNC_export$Approved_symbol)

# Load all data as Rds data
#raw_data = readRDS("./data/raw_data_paralog.Rds")
#Paraloc_data = readRDS("../../Paraloc_data_paraloc.Rds")


# generate_test_data ----

generate_test_data <- function() {
  
  # generate random input data for testing
  input<-data.frame(chr="1",pos="115256528",ref="T")
  input1<-data.frame(chr="1",pos="115256528",ref="T",alt="C")
  input2<-data.frame(chr="1",pos="115256528DD",ref="T",alt="G")
  input3<-data.frame(chr="3",pos="38592567",ref="TA",alt="B")
  input4<-data.frame(chr="21",pos="44592214",ref="WC",alt="T")
  input5<-data.frame(chr="21B",pos="47421902",ref="G",alt="A")
  input6<-data.frame(chr="XB",pos="70443591",ref="G",alt="A")

  
  input <- rbind(input1,input2,input3,input4,input5,input6)
  
  #var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  #var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
  var=var[nzchar(x=var)]
  input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
  input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2)
  
  #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
  
  
  return(input_data)
}


# function to test/validate variant input
validate_input <- function(input_data) {
  
  input_data <- suppressWarnings(tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F))
  
  #chr = c(as.character(1:22),'x','X','y','Y')
  input_data <- input_data[ ( input_data$CHR.query %in% c(as.character(1:22),'x','X','y','Y')  &  grepl("^\\d", suppressWarnings(as.numeric(input_data$POS.query))) & input_data$REF.query %in% c("A","C", "T", "G") & input_data$ALT.query %in% c("A","C", "T", "G") ) , ]
  
}


generate_test_data_2 <- function(var=NULL) {
  
  if (var==1){
  input<-data.frame(chr="1",pos="115256528",ref="T",alt="G")
  
  } else if (var==2){
    #3-38592567-T-A
    input<-data.frame(chr="3",pos="38592567",ref="T",alt="A")
    
    
  } else if (var==3){
    #21-44592214-C-T
    input<-data.frame(chr="21",pos="44592214",ref="C",alt="T")
    
  } else if (var==4){
    #X-70443591-G-A
    input<-data.frame(chr="X",pos="70443591",ref="G",alt="A")

  } else if (var==5){
    
    input1<-data.frame(chr="1",pos="115256528",ref="T",alt="G")
    input2<-data.frame(chr="3",pos="38592567",ref="T",alt="A")
    input3<-data.frame(chr="21",pos="44592214",ref="C",alt="T")
    input4<-data.frame(chr="X",pos="70443591",ref="G",alt="A")

    input <- rbind(input1,input2,input3,input4)
  }
  

  #var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  #var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
  var=var[nzchar(x=var)]
  input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
  input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2)
  
  #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
  
  
  return(input_data)
}

# validate_input ----

# function to test/validate variant input
validate_input <- function(input_data) {
  
  input_data <- suppressWarnings(tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F))
  
  #chr = c(as.character(1:22),'x','X','y','Y')
  input_data <- input_data[ ( input_data$CHR.query %in% c(as.character(1:22),'x','X','y','Y')  &  grepl("^\\d", suppressWarnings(as.numeric(input_data$POS.query))) & input_data$REF.query %in% c("A","C", "T", "G") & input_data$ALT.query %in% c("A","C", "T", "G") ) , ]
  
}


# edit_paralog_colnames ----

paraloc_colnames<- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  "var.query",
  "Gene.query",
  "Positions.paralog")

paralog_colnames <- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  "ALT.query",
  "var.query",
  "ID.query",
  "ClinVar.query",
  "Gene.query",
  "ENSG.query",
  "ENST.query",
  "cDNA.query",
  "Protein.query",
  "cDNA_position.query",
  "Protein_position.query",
  "AA.query",
  "Codons.query",
  "Para_Z_score.query",
  "CHR.paralog",
  "POS.paralog",
  "REF.paralog",
  "ALT.paralog",
  "var.paralog",
  "ID.paralog",
  "ClinVar.paralog",
  "Gene.paralog",
  "ENSG.paralog",
  "ENST.paralog",
  "cDNA.paralog",
  "Protein.paralog",
  "cDNA_position.paralog",
  "Protein_position.paralog",
  "AA.paralog",
  "Codons.paralog",
  "Para_Z_score.paralog")

homolog_colnames <- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  "ALT.query",
  "var.query",
  "ID.query",
  "ClinVar.query",
  "Gene.query",
  "ENSG.query",
  "ENST.query",
  "cDNA.query",
  "Protein.query",
  "cDNA_position.query",
  "Protein_position.query",
  "AA.query",
  "Codons.query",
  "Pfam_domain.query",
  "Pfam_pos.query",
  "CHR.homolog",
  "POS.homolog",
  "REF.homolog",
  "ALT.homolog",
  "var.homolog",
  "ID.homolog",
  "ClinVar.homolog",
  "Gene.homolog",
  "ENSG.homolog",
  "ENST.homolog",
  "cDNA.homolog",
  "Protein.homolog",
  "cDNA_position.homolog",
  "Protein_position.homolog",
  "AA.homolog",
  "Codons.homolog")


# lookup_vars ----

lookup_paralog <- function(input_data){
  
  paralog_out <- NULL
  
  for (i in 1:nrow(input_data)) {
    
    #i=1
    # check if ALT.query exists, if not skip and show only paraloc 
    #if (!is.na(input_data$ALT.query[i])) {
      query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
      #CMD_paralog<- paste0("tabix ", paralog_data, " ", query)
      #tabix_paralog <- system(command = paste0("tabix ", paralog_data, " ", query), intern = T,wait = T)
  
      #tabix_paralog <- system(command = paste0("tabix data/paralog_data_sorted.txt.gz ", query), intern = T,wait = T)
      tabix_paralog_extra <- suppressMessages(suppressWarnings(system(command = paste0("tabix data/paralog_data.txt.gz ", query), intern = T,wait = T)))

      
      #pg1 <- separate(as.data.frame(unlist(tabix_paralog)),1,sep = "\t", into =  paralog_colnames)
      #pg1 <- pg1[(pg1$REF.query==input_data$REF.query[i] & pg1$ALT.query == input_data$ALT.query[i]),]
      
      #pg1 <- as.data.frame(stringr::str_split_fixed(tabix_paralog, pattern = "\t", length(paralog_colnames)), stringsAsFactors = F)
      pg1 <- as.data.frame(stringr::str_split_fixed(tabix_paralog_extra, pattern = "\t", length(paralog_colnames)), stringsAsFactors = F)
      
      #colnames(pg1) <- paralog_colnames
      
      
      if (is.null(paralog_out)){
        paralog_out = pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),]
      } else {
        paralog_out = rbind(paralog_out, pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),])
      }
    #}
  }
  #paralog_out = rbind(paralog_out, pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),])
  
  #colnames(paralog_out) <- paralog_colnames
  colnames(paralog_out) <- paralog_colnames
  
  # sort df to avoid sorting later
  #chrOrder<-c(1:22,"X","Y")
  paralog_out <- paralog_out[order(factor(paralog_out$CHR.query , levels = c(1:22,"X","Y")), 
                                   as.numeric(paralog_out$POS.query),
                                   factor(paralog_out$CHR.paralog , levels = c(1:22,"X","Y")), 
                                   as.numeric(paralog_out$POS.paralog)), ]
  

  return(paralog_out)
}


lookup_paraloc <- function(input_data){
  
  paraloc_out <- NULL
  
  for (i in 1:nrow(input_data)) {
    
    query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
    #CMD_paraloc<- paste0("tabix ", paraloc_data, " ", query)
    #tabix_paraloc <- system(command =  paste0("tabix ", paraloc_data, " ", query), intern = T, wait = T)
    
    # tabix_paraloc <- system(command =  paste0("tabix data/paraloc_data_sorted.txt.gz " , query), intern = T, wait = T)
    tabix_paraloc <-suppressMessages(suppressWarnings(system(command =  paste0("tabix data/paraloc_chr/paraloc_data_chr", input_data$CHR.query[i] ,".txt.gz " , query), intern = T, wait = T)))
    
    
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

  return(unique(paraloc_out))
  
}




lookup_homolog <- function(input_data){
  
  homolog_out <- NULL
  #input_data <- validate_input(generate_test_data_2(var = 1))

  for (i in 1:nrow(input_data)) {
    
    #i=1
    query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
    tabix_homolog <- suppressMessages(suppressWarnings(system(command = paste0("tabix data/homolog_data.txt.gz ", query), intern = T,wait = T)))
    #homo <- system(command = paste0("tabix data/homolog_data.txt.gz ", query), intern = T,wait = T)
    #pg1 <- separate(as.data.frame(unlist(tabix_paralog)),1,sep = "\t", into =  paralog_colnames)
    #pg1 <- pg1[(pg1$REF.query==input_data$REF.query[i] & pg1$ALT.query == input_data$ALT.query[i]),]
    hg1 <- as.data.frame(stringr::str_split_fixed(tabix_homolog, pattern = "\t", length(homolog_colnames)), stringsAsFactors = F)

    if (is.null(homolog_out)){
      homolog_out = hg1[(hg1$V3==input_data$REF.query[i] & hg1$V4 == input_data$ALT.query[i]),]
    } else {
      homolog_out = rbind(homolog_out, hg1[(hg1$V3==input_data$REF.query[i] & hg1$V4 == input_data$ALT.query[i]),])
    }
  }
  
  colnames(homolog_out) <- homolog_colnames
  
  # sort df to avoid sorting later
  #chrOrder<-c(1:22,"X","Y")
  homolog_out <- homolog_out[order(factor(homolog_out$CHR.query , levels = c(1:22,"X","Y")), 
                                   as.numeric(homolog_out$POS.query),
                                   factor(homolog_out$CHR.homolog , levels = c(1:22,"X","Y")), 
                                   as.numeric(homolog_out$POS.homolog)), ]
  
  
  return(homolog_out)
}




# predict_output ----



predict_output_tabix = function(input_data){


  #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F) 
  # input_data <- validate_input(generate_test_data_2(1))
  
    paralog <- lookup_paralog(input_data)
    paraloc <- lookup_paraloc(input_data)
    homolog <- lookup_homolog(input_data)
    

  return(list("paralog" = paralog, "paraloc" = paraloc, "homolog" = homolog))

}  


# result <- predict_output_tabix(input_data = input_line)

  
predict_output = function(input_data){
  #### select through RDS objects ####
  paralog = raw_data[raw_data$var %in%  input_data$mutation,]
  paraloc = Paraloc_data[Paraloc_data$var %in% input_data$paraloc,]
  
  # will change subset to dplyr::select to change colnames at the same time
  paralog= dplyr::select(paralog_output,
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
  paralog = dplyr::mutate_if(paralog, is.numeric, as.character)
  
  return(list("paralog" = paralog, "paraloc" = paraloc))

}



# check_upload_file ----


# function to check if uploaded variants file is valid
check_upload_file = function(inFile) {
  # very hacky way to read in vcf
  # read 1st line only to get number of cols to check if its txt or vcf format
  # when using colClasses in read.table the columns set to NULL are completely ignored
  # input_1row = ncol(read.table(inFile$datapath,nrows = 1 ))
  
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

# add_URLs ---- 

# example data
# input_data_2 <- generate_test_data_2(1)
# # input_data_2 <- rbind(generate_test_data_2(1) ,generate_test_data_2(2))
# result_paraloc <- predict_output_tabix(validate_input(input_data_2))$paraloc


# function to tap into Ensembl API to retreive transcriptid, gene id and protein position for paraloc position varint.
query_paraloc_API <- function(para_split_df_long) {
  
  
  for (i in 1:nrow(para_split_df_long)) {
  
    # paraloc_gene <- para_split_df_long$Gene.paraloc[2]
    # para_split_df_long <- para_split_df_long[2,]
    # url_string <- paste0("https://grch37.rest.ensembl.org/vep/human/region/",para_split_df_long$chr.paraloc,":",para_split_df_long$AA_pos.paraloc,":1/G?canonical=1")

    url_string <- paste0("https://grch37.rest.ensembl.org/vep/human/region/",para_split_df_long$chr.paraloc[i],":",para_split_df_long$AA_pos.paraloc[i],":1/G?canonical=1")
    
    API_out <- content(GET(url_string , content_type("text/plain")))[[1]] #[[1]]$transcript_consequences
    
    #API_out_1 <- API_out[[1]]
    #transcript_consequences <- API_out$transcript_consequences
    
    for ( l in c(1:length(API_out[["transcript_consequences"]]))) {
      
      #transcript_consequences <-  transcript_consequences[[l]]
      
      if ((!is.null(API_out[["transcript_consequences"]][[l]][["canonical"]]) && API_out[["transcript_consequences"]][[l]][["canonical"]]==1) && ((!is.null(API_out[["transcript_consequences"]][[l]][["gene_symbol"]]) && API_out[["transcript_consequences"]][[l]][["gene_symbol"]]==para_split_df_long$Gene.paraloc[i]))){
        para_split_df_long$ENST.paraloc[i] <- API_out[["transcript_consequences"]][[l]]["transcript_id"]
        para_split_df_long$ENSG.paraloc[i] <- API_out[["transcript_consequences"]][[l]]["gene_id"]
        para_split_df_long$Protein_position.paraloc[i] <- API_out[["transcript_consequences"]][[l]]["protein_start"]
      }
      
      # if ((!is.null(transcript_consequences[[l]]$canonical) && transcript_consequences[[l]]$canonical==1) && ((!is.null(transcript_consequences[[l]]$gene_symbol) && transcript_consequences[[l]]$gene_symbol==para_split_df_long$Gene.paraloc[i]))){
      #   para_split_df_long$ENST.paraloc[i] <- transcript_consequences[[l]]$transcript_id
      #   para_split_df_long$ENSG.paraloc[i] <- transcript_consequences[[l]]$gene_id
      #   para_split_df_long$Protein_position.paraloc[i] <- transcript_consequences[[l]]$protein_start
      # }
      # 
    }
    
    
  }
  return(para_split_df_long)
  
}




# function to add ensembl URLs to paralg positions in genes
add_paraloc_URL_new = function(result_paraloc) {
  
  if (nrow(result_paraloc)!=0){
    
    #embl_api <- "https://grch37.rest.ensembl.org"
    para_split_df <- NULL
    
    
    result_paraloc$Positions.paralog <- str_split(result_paraloc$Positions.paralog , pattern = ", ", simplify = F)
    
    for (i in 1:nrow(result_paraloc)){ 
      # get all positions from list and apply add_URLs function to every position
      # then paste back as string
      
      # uncomment to use inside the loop
      # para_split_df_long <- suppressMessages(suppressWarnings(cbind(result_paraloc[1,1:5] ,
      #                            separate(data = (data.frame(result_paraloc$Positions.paralog[1]) %>% rename(para_split=1)),col = "para_split" ,into = c("Gene.paraloc", "chr.paraloc", "AA_pos.paraloc","AA.paraloc"),sep = " " ))))


      para_split_df_long <- suppressMessages(suppressWarnings(cbind(result_paraloc[i,1:5] , 
                            separate(data = (data.frame(result_paraloc$Positions.paralog[i]) %>% rename(para_split=1)),col = "para_split" ,into = c("Gene.paraloc", "chr.paraloc", "AA_pos.paraloc","AA.paraloc"),sep = " " ))))
      
      para_split_df_long$ENST.paraloc <- NA
      para_split_df_long$ENSG.paraloc <- NA
      para_split_df_long$Protein_position.paraloc <- NA
      
      # run the query to Ensembl API to get transcript, gene and protein position for paraloc
      ###
      # uncomment bellow line to query API
      #para_split_df_long_API <- query_paraloc_API(para_split_df_long)  ############## revert API call 
      para_split_df_long_API <- para_split_df_long
      

      # example code to tap into Ensembl API
      #para_split_df_long$Codons <- content(GET(paste(embl_api, ext[i], sep = ""), content_type("text/plain")))
    
      
      if (is.null(para_split_df)){
        para_split_df = para_split_df_long_API
      } else {
        para_split_df = rbind(para_split_df, para_split_df_long_API)
      }
      
    }

      
    #Ensembl Gene.query for paraloc
    para_split_df$Gene.query<- ifelse(!is.na(para_split_df$Gene.query),
                                      (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(para_split_df$Gene.query)]), "' target='_blank'>", para_split_df$Gene.query, "</a>")),
                                      "-")

    #Ensembl gene for paraloc
    para_split_df$Gene.paraloc<- ifelse(!is.na(para_split_df$Gene.paraloc),
                                        (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(para_split_df$Gene.paraloc)]), "' target='_blank'>", para_split_df$Gene.paraloc, "</a>")),
                                        "-")

    #Ensembl ENST.paraloc
    para_split_df$ENST.paraloc<- ifelse(!is.na(para_split_df$ENST.paraloc),
                                        (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",para_split_df$ENSG.paraloc,";t=",para_split_df$ENST.paraloc), "' target='_blank'>", para_split_df$ENST.paraloc, "</a>")),
                                        "-")

    #para_split_df <- para_split_df[1:8,10,12,9]
    para_split_df <- para_split_df[order(factor(para_split_df$CHR.query , levels = c(1:22,"X","Y")),
                                         as.numeric(para_split_df$POS.query),
                                         factor(para_split_df$chr.paraloc , levels = c(1:22,"X","Y")) #,as.numeric(para_split_df$AA_pos.paraloc) # its a character so cannot be ordered properly
                                         ), c(1:8,10,12,9) ]
    #print(names(para_split_df))
    para_split_df <- para_split_df[,c(1:5,11,6:10)]
  } else { 
    para_split_df <- NULL
  }
  
  return(para_split_df)
}




# function to add ensembl URLs to paralg positions in genes
add_paralog_URL = function(result) {
  

    #ClinVarID paralog URL
    result$ID.paralog<- ifelse(#!is.na(result$ID.paralog),
      result$ID.paralog!="NA",
      (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.paralog,"/"), "' target='_blank'>", result$ID.paralog, "</a>")),
      "-")
    
    #ClinVarID query URL
    result$ID.query<- ifelse(#!is.na(result$ID.query),
      result$ID.query!="NA",
      (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.query,"/"), "' target='_blank'>", result$ID.query, "</a>")),
      "-")
    #print(paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene.query)],";g1=",map[unlist(result$SYMBOL.paralog)]))
    #Ensembl alignment URL
    # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
    result$Ensembl_alignment_link<- ifelse(!is.na(result$ENSG.paralog), 
                                           (paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",result$ENSG.query,";g1=",result$ENSG.paralog), "' class='btn btn-default btn-sm btn-block active' target='_blank'>alignment</a>")) , 
                                           "-") 
    
    # https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000213281;t=ENST00000369535
    #Ensembl ENST.query
    result$ENST.query<- ifelse(!is.na(result$ENST.query), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",result$ENSG.query,";t=",result$ENST.query), "' target='_blank'>", result$ENST.query, "</a>")),
                               "-")
    
    #Ensembl ENST.paralog
    result$ENST.paralog<- ifelse(!is.na(result$ENST.paralog), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",result$ENSG.paralog,";t=",result$ENST.paralog), "' target='_blank'>", result$ENST.paralog, "</a>")),
                               "-")
    
    #HGNC Gene.query
    result$Gene.query<- ifelse(!is.na(result$Gene.query), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.query), "' target='_blank'>", result$Gene.query, "</a>")),
                               "-")
    
    #HGNC Gene.paralog
    result$Gene.paralog<- ifelse(!is.na(result$Gene.paralog), 
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.paralog), "' target='_blank'>", result$Gene.paralog, "</a>")),
                                 "-")
    
    result <- cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', result )
    
    return(result)
    
} 

# function to add ensembl URLs to homolog positions in genes
add_homolog_URL = function(result) {
  
  
  result$Pfam_pos.query <- strsplit(result$Pfam_pos.query,split = "=")[[1]][2]
  
  #ClinVarID paralog URL
  result$ID.homolog<- ifelse(#!is.na(result$ID.paralog),
    result$ID.homolog!="NA",
    (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.homolog,"/"), "' target='_blank'>", result$ID.homolog, "</a>")),
    "-")
  
  #ClinVarID query URL
  result$ID.query<- ifelse(#!is.na(result$ID.query),
    result$ID.query!="NA",
    (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",result$ID.query,"/"), "' target='_blank'>", result$ID.query, "</a>")),
    "-")
  #print(paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(result$Gene.query)],";g1=",map[unlist(result$SYMBOL.paralog)]))
  #Ensembl alignment URL
  # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
  # result$Ensembl_alignment_link<- ifelse(!is.na(result$ENSG.paralog), 
  #                                        (paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",result$ENSG.query,";g1=",result$ENSG.paralog), "' class='btn btn-default btn-sm btn-block active' target='_blank'>alignment</a>")) , 
  #                                        "-") 
  
  # https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000213281;t=ENST00000369535
  #Ensembl ENST.query
  result$ENST.query<- ifelse(!is.na(result$ENST.query), 
                             (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",result$ENSG.query,";t=",result$ENST.query), "' target='_blank'>", result$ENST.query, "</a>")),
                             "-")
  
  #Ensembl ENST.paralog
  result$ENST.homolog<- ifelse(!is.na(result$ENST.homolog), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",result$ENSG.homolog,";t=",result$ENST.homolog), "' target='_blank'>", result$ENST.homolog, "</a>")),
                               "-")
  
  #HGNC Gene.query
  result$Gene.query<- ifelse(!is.na(result$Gene.query), 
                             (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.query), "' target='_blank'>", result$Gene.query, "</a>")),
                             "-")
  
  #HGNC Gene.paralog
  result$Gene.homolog<- ifelse(!is.na(result$Gene.homolog), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.homolog), "' target='_blank'>", result$Gene.homolog, "</a>")),
                               "-")
  
  result <- cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', result )
  return(result)
  
} 


# DT_colnames ----

paralog_DT_colnames <- c('Chr.query' = 'CHR.query',
                         'Pos.query' = 'POS.query',
                         'REF.query' = 'REF.query',
                         'ALT.query' = 'ALT.query',
                         'Query variant' = 'var.query',
                         'ClinVar.query' = 'ID.query',
                         'ClinVar Class.query' = 'ClinVar.query',
                         'Gene.query' = 'Gene.query',
                         'ENST.query' = 'ENST.query',
                         'ENSG.query' = 'ENSG.query',
                         'cDNA.query' = 'cDNA.query',
                         'Protein.query' = 'Protein.query',
                         'cDNA position.query' = 'cDNA_position.query',
                         'Protein position.query' = 'Protein_position.query',
                         'AA.query' = 'AA.query',
                         'Codons.query' = 'Codons.query',
                         'Para_Z Score.query'='Para_Z_score.query',
                         'Chr' = 'CHR.paralog',
                         'Pos' = 'POS.paralog',
                         'REF' = 'REF.paralog',
                         'ALT' = 'ALT.paralog',
                         'Paralogous variant' = 'var.paralog',
                         'ClinVar ID' = 'ID.paralog',
                         'ClinVar Class' = 'ClinVar.paralog',
                         'Gene' = 'Gene.paralog',
                         'ENSG' = 'ENSG.paralog',
                         'ENST' = 'ENST.paralog',
                         'cDNA' = 'cDNA.paralog',
                         'Protein' = 'Protein.paralog',
                         'cDNA position' = 'cDNA_position.paralog',
                         'Protein position' = 'Protein_position.paralog',
                         'AA' = 'AA.paralog',
                         'Codons' = 'Codons.paralog',
                         'Para_Z Score' = 'Para_Z_score.paralog',
                         'Ensembl alignment' = 'Ensembl_alignment_link')



paraloc_DT_colnames <- c(
  'Query variant' = 'var.query',
  'Query gene' = 'Gene.query',
  'Query Residue' = 'AA.paraloc',
  'Paralogous gene' = 'Gene.paraloc',
  'Chromosome' = 'chr.paraloc',
  'AA Position' = 'AA_pos.paraloc',
  #'AA Residue' = 'AA.paraloc',
  'ENST' = 'ENST.paraloc',
  #'ENSG' = 'ENSG.paraloc',
  'Protein positions' = 'Protein_position.paraloc'
  )

homolog_DT_colnames <- c('Chr.query' = 'CHR.query',
                         'Pos.query' = 'POS.query',
                         'REF.query' = 'REF.query',
                         'ALT.query' = 'ALT.query',
                         'Query variant' = 'var.query',
                         'ClinVar.query' = 'ID.query',
                         'ClinVar Class.query' = 'ClinVar.query',
                         'Gene.query' = 'Gene.query',
                         'ENST.query' = 'ENST.query',
                         'ENSG.query' = 'ENSG.query',
                         'cDNA.query' = 'cDNA.query',
                         'Protein.query' = 'Protein.query',
                         'cDNA position.query' = 'cDNA_position.query',
                         'Protein position.query' = 'Protein_position.query',
                         'AA.query' = 'AA.query',
                         'Codons.query' = 'Codons.query',
                         'Pfam domain' = 'Pfam_domain.query',
                         'Pfam position' = 'Pfam_pos.query',
                         'Chr' = 'CHR.homolog',
                         'Pos' = 'POS.homolog',
                         'REF' = 'REF.homolog',
                         'ALT' = 'ALT.homolog',
                         'Homologous variant' = 'var.homolog',
                         'ClinVar ID' = 'ID.homolog',
                         'ClinVar Class' = 'ClinVar.homolog',
                         'Gene' = 'Gene.homolog',
                         'ENSG' = 'ENSG.homolog',
                         'ENST' = 'ENST.homolog',
                         'cDNA' = 'cDNA.homolog',
                         'Protein' = 'Protein.homolog',
                         'cDNA position' = 'cDNA_position.homolog',
                         'Protein position' = 'Protein_position.homolog',
                         'AA' = 'AA.homolog',
                         'Codons' = 'Codons.homolog')



# childrow_JS_callback ----

childrow_JS_callback_paralog <- c("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<table style=\"border-spacing:50px\" style=\"width:50%\" cellpadding=\"50px\" style=\"padding-left:50px\">'+
      '<thead>'+
    //  '<tr>'+
    //        '<th colspan=\"4\">Query variant</td>'+
    //  '</tr>'+
    //    '<tr>'+
    //        '<th>Query variant</th>'+
    //        '<th>ClinVar ID</th>'+
    //        '<th>ClinVar Class</th>'+
    //        '<th>Gene</th>'+
    //        '<th>HGVS</th>'+
    //        '<th>Codons</th>'+
    //        '<th>Para_Z Score</th>'+
    //    '</tr>'+
    //   '</thead>'+
    //   '<tbody>'+
    //    '<tr>'+
    //        '<td>'+ d[1] + '-' + d[2] + '-' + d[3] + '-' + d[4] +'</td>'+
    //        '<td>'+ d[6] +'</td>'+
    //        '<td>'+ d[7] +'</td>'+
    //        '<td>'+ d[8] +'</td>'+
    //        '<td>'+ d[10] + '(' + d[9] + '):' + d[11] + ' (' + d[12] + ')' +'</td>'+
    //        '<td>'+ d[16] +'</td>'+
    //        '<td>'+ d[17] +'</td>'+
    //    '</tr>'+
    //  '</tbody>'+
    //  '<tr>'+
    //      '<th>Query variant</th>'+
    //      '<td>'+ d[1] + '-' + d[2] + '-' + d[3] + '-' + d[4] +'</td>'+
    //  '</tr>'+
        '<tr>'+
            '<th>ClinVar ID</th>'+
            '<td>'+ d[6] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>ClinVar class</th>'+
            '<td>'+ d[7] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Gene</th>'+
            '<td>'+ d[8] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>HGVS</th>'+
    //      '<td>'+ d[9] + '(' + d[8] + '):' + d[10] + ' (' + d[11] + ')' +'</td>'+
            '<td>'+ d[10] + ':' + d[11] + ' (' + d[12] + ')' +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Codons</th>'+
            '<td>'+ d[16] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Para_Z Score</th>'+
            '<td>'+ d[17] +'</td>'+
        '</tr>'+
    '</table>';
  };
  table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
  //    td.html('<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css\"> <i class=\"fa fa-plus-square fa-lg\"></i>');
    } else {
      row.child(format(row.data())).show();
  //    td.html('<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css\"> <i class=\"fa fa-minus-square fa-lg\"></i>');
    }
  }
);")


childrow_JS_callback_homolog <- c("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<table style=\"border-spacing:50px\" style=\"width:50%\" cellpadding=\"50px\" style=\"padding-left:50px\">'+
      '<thead>'+
    //  '<tr>'+
    //      '<th>Query variant</th>'+
    //      '<td>'+ d[1] + '-' + d[2] + '-' + d[3] + '-' + d[4] +'</td>'+
    //  '</tr>'+
        '<tr>'+
            '<th>ClinVar ID</th>'+
            '<td>'+ d[6] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>ClinVar class</th>'+
            '<td>'+ d[7] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Gene</th>'+
            '<td>'+ d[8] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>HGVS</th>'+
    //      '<td>'+ d[9] + '(' + d[8] + '):' + d[10] + ' (' + d[11] + ')' +'</td>'+
            '<td>'+ d[10] + ':' + d[11] + ' (' + d[12] + ')' +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Codons</th>'+
            '<td>'+ d[16] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Pfam domain</th>'+
            '<td>'+ d[17] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Pfam domain position</th>'+
            '<td>'+ d[18] +'</td>'+
        '</tr>'+
    '</table>';
  };
  table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
  //    td.html('<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css\"> <i class=\"fa fa-plus-square fa-lg\"></i>');
    } else {
      row.child(format(row.data())).show();
  //    td.html('<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css\"> <i class=\"fa fa-minus-square fa-lg\"></i>');
    }
  }
);")


####### draw proteins example ######

# "P51787 O43526 O43525 P56696 Q9NR82"
# "Q04206 Q01201 Q04864 P19838 Q00653"


####################################

# example data
# input_data_2 <- generate_test_data_2(4) 
#input_data_2 <- rbind(generate_test_data_2(1) ,generate_test_data_2(4))
# 
#result <- predict_output_tabix(validate_input(input_data_2))$paralog
#result <- add_paraloc_URL_new(predict_output_tabix(validate_input(input_data_2))$paraloc)

# prepare draw protein dataframe ----


# function to order protein accession numbers to dispaly the query prot first
query_protein_API <- function(query_proteins) {
  
  # example query_proteins
  #query_proteins <- "P51787"
  #query_proteins <- "P51787 O43526 O43525 P56696 Q9NR82"
  
  # GET protein info from UniProt API
  #prot_data <- feature_to_dataframe(get_features(query_proteins))
  
  
  # get pfam info from Pfam API
  # ###changed the function get_Pfam to include order col
  pfam_data <- get_Pfam(query_proteins)
  
  # left join UniProt order 
  #pfam_data2 <- pfam_data %>% left_join(prot_data %>% select(accession, taxid, order) %>% unique.array(),by = "accession")
  pfam_data$begin <- as.numeric(pfam_data$begin)
  
  
  # join Uniprot and Pfam data
  #prot_data <- rbind(prot_data, pfam_data)
  
  #return(prot_data)
  return(pfam_data)
  
}



# 
get_prot_data <- function(result) {
  

  result$UniProt_ID.paralog <- map_HGNC[unlist(result$Gene.paralog)]
  result$UniProt_ID.query <- map_HGNC[unlist(result$Gene.query)]
  
  #query_proteins <- paste(map_HGNC[unlist(unique(result$Gene.query))], paste(unique(result$UniProt_ID.paralog),collapse = " "))
  query_proteins_query  <- paste(map_HGNC[unlist(unique(result$Gene.query))])
  query_proteins_paralog<- paste(unique(result$UniProt_ID.paralog),collapse = " ")
  
  
  
  
  query_prot <- query_protein_API(query_proteins = query_proteins_query)
  paralog_prot <-  query_protein_API(query_proteins = query_proteins_paralog)
  
  # reorder query prot so it comes up as first in the graph
  query_prot$order <- as.numeric(max(paralog_prot$order)+1)
  
  #rbind query and paralog prot
  prot_data <- rbind(query_prot,paralog_prot)
  
  # left join protein data
  prot_data <- prot_data %>% left_join(rbind(result %>% select("UniProt_ID" = UniProt_ID.query, "Gene" = Gene.query, "ENSG" = ENSG.query, "Protein_position" = Protein_position.query) %>% unique.array(),
                                             result %>% select("UniProt_ID" = UniProt_ID.paralog, "Gene" = Gene.paralog, "ENSG" = ENSG.paralog, "Protein_position" = Protein_position.paralog) %>% unique.array()),
                                       by = c("accession"= "UniProt_ID"))


  
  
  
  prot_data$Protein_position <- as.numeric(prot_data$Protein_position)
  #query_HGVS
  prot_data$HGVS.query <- paste(unique(result$ENST.query),"(", unique(result$Gene.query),"):",unique(result$cDNA.query), " (",unique(result$Protein.query),")", sep = "")
  #query_var_ID
  prot_data$query_var_ID <- paste(unique(result$var.query))

  #Add Gene URL
  prot_data$Gene<- ifelse(
    prot_data$Gene!="NA",
    (paste0("<a href='", paste0("http://pfam.xfam.org/protein/",prot_data$accession,"/"), "' target='_blank'>", prot_data$Gene, "</a>")),
    prot_data$Gene)
  
  
  prot_data <- separate(data = prot_data,col = description, into = "description", extra = "drop", sep = ";")
  
  return(prot_data)
}
#####


# draw the prot graph ----
draw_prot_data_plotly <- function(input_data) {
  # get prot data
  
  # #test
  # input_data_2 <- generate_test_data_2(4)
  # input_data_2 <- rbind(generate_test_data_2(1) ,generate_test_data_2(4))
  # # 
  prot_data <- list()
  for ( i in unique(input_data$var.query)) {
    
    prot_slice <- input_data[input_data$var.query == i,]
    
    
    name <- paste(i)
    tmp <- get_prot_data( prot_slice)
    prot_data[[name]] <- tmp
    
  }

  
  # subplot(lapply(prot_data, draw_plotly_graph),margin = 0.05,nrows = length(prot_data),)
  #fig <- lapply(prot_data, draw_plotly_graph)
  fig <- lapply(prot_data, draw_plotly_graph_PFAM)
  
  
  
  # #test
  #
  # input_data_2 <- generate_test_data_2(3)
  # input_data <- predict_output_tabix(validate_input(input_data_2))$paralog
  # prot_data <- get_prot_data(input_data[input_data$var.query == unique(input_data$var.query)[1],])# %>% dplyr::arrange(desc(order))
  #
  #  fig <- draw_plotly_graph(prot_data = prot_data)
  # fig
  # #
  
}



# draw plotly graph PFAM positions ----

draw_plotly_graph_PFAM <- function(prot_data, showlegend = T) {
  

  fig <- plot_ly(prot_data)
  
  # Chains
  fig <- add_bars(fig, data = prot_data[prot_data$type == "PFAM_name",],x = ~c((begin)-(end)), y = ~Gene, base = ~(end-Protein_position),
                  width = 0.02, orientation = 'h', showlegend = F ,name = ~Gene, marker = list(color = toRGB("gray50")), 
                  hoverinfo = "name+text") # color = 'rgba(50, 171, 96, 0.6)'
  # Pfam domains
  fig <- add_bars(fig, data = prot_data[prot_data$type == "PFAM",], x = ~c((begin)-(end)), y = ~Gene, base = ~(end-Protein_position) ,#type = 'bar',
                   width = 0.4, orientation = 'h', showlegend = T, name = ~description, 
                   hoverinfo = "name+text")
  # HGVS
  fig <- add_annotations(fig, data = prot_data,x = 0,y = 1, yref = "paper", yanchor = "bottom", showarrow = F, font = list(size = 14),
                         text = paste0(unique(prot_data$HGVS.query)) )  
  
  # legend , title and margins
  fig <- layout(fig,
                barmode = 'overlay', showlegend = T,
                title = list(text = paste0(unique(prot_data$query_var_ID)),font = list(size = 16), 
                             xref = "paper",yref = "paper",xanchor = "left", x = 0 ),
                legend = list(itemdoubleclick = "toggle", title = list(text = "Pfam Domains",font = list(size = 14))),
                yaxis = list(title = "",autorange = T, showgrid = F, showline = F, showticklabels = T) ,margin= c(0, 0.95),# will set the total height of the plot 
                xaxis = list(title = "",autorange = T, showgrid = T, showline = F, showticklabels = F,zeroline = T, zerolinewidth = 3)) 
  
  # plotly toolbox
  fig <- plotly::config(fig, displayModeBar = T, 
                        modeBarButtonsToRemove = list("pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d","resetScale2d","hoverClosestCartesian", "hoverCompareCartesian", "hoverClosestGl2d", "toggleSpikelines"), 
                        toImageButtonOptions = list(format = "png",width = "1800", height = "800"), displaylogo = F )
  
  # fig

  
  return(fig)
  

}
#####


# draw plotly graph ----

draw_plotly_graph <- function(prot_data, showlegend = T, title ) {
  
  fig <- plot_ly(prot_data)
  
  # Chains
  if ("CHAIN" %in% prot_data$type ) {
    fig <- add_bars(fig, data = prot_data[prot_data$type == "CHAIN",],x = ~c((begin)-(end)), y = ~reorder(Gene, order), base = ~(end-Protein_position),
                    width = 0.02, orientation = 'h', showlegend = F , legendgroup = "Chains", # yaxis = ~Gene,  #xaxis = Gene
                    marker = list(color = toRGB("gray50")), name = ~Gene, hoverinfo = "name+text") # color = 'rgba(50, 171, 96, 0.6)'
  }
  
  
  # Folding
  if ("HELIX" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "HELIX",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.2, orientation = 'h', showlegend = F, name = ~type, legendgroup = "Folding",
                     opacity = 0.1 , hoverinfo = "name+text")
  }
  
  if ("STRAND" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "STRAND",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.2, orientation = 'h', showlegend = F, name = ~type,  legendgroup = "Folding",
                     opacity = 0.1 , hoverinfo = "name+text")
  }
  
  if ("TURN" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "TURN",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.2, orientation = 'h', showlegend = F, name = ~type,  legendgroup = "Folding",
                     opacity = 0.1 , hoverinfo = "name+text")
  }
  
  # Repeat
  if ("REPEAT" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "REPEAT",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.1, orientation = 'h', showlegend = F, name = ~description,  legendgroup = 'Repeat',
                     opacity = 0.2, marker = list(color = toRGB("gray50"),
                                                  line = list(color = toRGB("gray20"), width = 2)),
                     hoverinfo = "name+text")
  }
  
  # Region
  if ("REGION" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "REGION",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  # Domain
  if ("DOMAIN" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "DOMAIN",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  
  # Topo_Domain
  if ( ("TOPO_DOM" %in% prot_data$type ) | ("TRANSMEM" %in% prot_data$type)) {
    fig <- add_trace(fig, data = prot_data[(prot_data$type == "TOPO_DOM" | prot_data$type == "TRANSMEM") ,], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  
  
  # Motif
  if ("MOTIF" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "MOTIF",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  
  
  # PFam
  if ("PFAM" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "PFAM",], x = ~c((begin)-(end)), y = ~reorder(Gene, order),type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'PFam', 
                     hoverinfo = "name+text") 
  }
  
  # # Phosphorilation NOT WORKING
  # if ("MOD_RES" %in% prot_data$type) {
  #   fig <- add_trace(fig, data = phospho_site_info(prot_data), x = ~c((begin)-(end)), y = ~Gene,type = 'scatter', base = ~(end-Protein_position),
  #                    width = 0.04, orientation = 'h', showlegend = T, name = "Phosphorilation",  # legendgroup = 'PFam', 
  #                    hoverinfo = "name+text") 
  # }
  
  fig <- add_annotations(fig, data = prot_data,x = 0,y = 1, yref = "paper", yanchor = "bottom", showarrow = F, font = list(size = 12),
                         text = paste0(unique(prot_data$HGVS.query)) )  
  
  fig <- layout(fig, 
                yaxis = list(title = "",autorange = T, showgrid = F, showline = F, showticklabels = T),   margin= c(0, 0.95), # will set the total height of the plot 
                title = list(text = paste0(unique(prot_data$query_var_ID)),
                              font = list(size = 16), xref = "paper",xref = "paper",xanchor = "left", x = 0 ),
                xaxis = list(title = "",autorange = T, showgrid = T, showline = F, showticklabels = F,zeroline = T, zerolinewidth = 3),
                barmode = 'overlay', showlegend = showlegend,legend = list(traceorder = "grouped", 
                                                                           itemdoubleclick = "toggle", 
                                                                           tracegroupgap = 20,
                                                                           title = list(text = "Grouped Domain Annotations",font = list(size = 14)))) 
  fig <- plotly::config(fig, displayModeBar = T, 
                modeBarButtonsToRemove = list("pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d","resetScale2d","hoverClosestCartesian", "hoverCompareCartesian", "hoverClosestGl2d", "toggleSpikelines"), 
                toImageButtonOptions = list(format = "png",width = "1800", height = "800"), displaylogo = F )
  
  
  return(fig)
}

######

# get pfam data ----
get_Pfam <- function (proteins_acc) 
{
  proteins_acc_url <- gsub(" ", "%2C", proteins_acc)
  
  baseurl <- "http://pfam.xfam.org/protein/api/"
  url <- paste0(baseurl, proteins_acc_url)
  
  pfam_eg <- suppressWarnings(suppressMessages(read_delim(url, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)))[,1:5] #%>% select()
  
  if (nrow(pfam_eg)==0) {
    print(paste("An error has occured while trying to connect to Pfam API."))
    #print("Pfam download has worked.")
  }
  # else {
  #   print(paste("An error has occured."))
  # }
  
  order=0
  pfam_eg[c("type", "description", "begin", "end", "length", "accession", "entryName", "order")] <- NA
  for (i in 1:nrow(pfam_eg)) {
    # set up cols
    if (pfam_eg[i,1]=="P") {
      pfam_eg$type[i] <- "PFAM_name"
      pfam_eg$description[i] <- pfam_eg$X4[i]
      pfam_eg$begin[i] <- "1"
      pfam_eg$end[i] <- pfam_eg$X5[i] 
      pfam_eg$length[i] <- pfam_eg$X5[i] -1
      pfam_eg$accession[i] <- pfam_eg$X2[i]
      pfam_eg$entryName[i] <- pfam_eg$X3[i]
      
      order <- order + 1
      pfam_eg$order[i] <- order
      
    } else if (pfam_eg[i,1]=="A") {
      pfam_eg$type[i] <- "PFAM"
      pfam_eg$description[i] <- pfam_eg$X2[i]
      pfam_eg$begin[i] <- pfam_eg$X4[i]
      pfam_eg$end[i] <- pfam_eg$X5[i]
      pfam_eg$length[i] <- pfam_eg$X5[i] - as.numeric(pfam_eg$X4[i])
      pfam_eg$entryName[i] <- pfam_eg$X3[i]
      
      
      pfam_eg$order[i] <- order 
      
    }
    if (is.na(pfam_eg$accession[i])){
      pfam_eg$accession[i] <- pfam_eg$accession[i-1]
    }
    
  }

  return(pfam_eg[,6:13])
}

######




