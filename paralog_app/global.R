library(shiny)
library(DT)
library(shinythemes)
library(tidyverse)
#library(microbenchmark)
#library(drawProteins)
library(plotly)
library(grid)




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




#paralog_index<- read.delim("data/paralog_data_index.txt", stringsAsFactors = F)
#saveRDS(paralog_index, "data/paralog_data_index.Rds")
#paralog_index <- readRDS("data/paralog_data_index.Rds")

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

generate_test_data_1 <- function() {
  
  # generate random input data for testing
  input<-data.frame(chr="1",pos="115256528",ref="T")

  #var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  var = paste(input$chr,input$pos,input$ref,sep = " ")
  #var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
  var=var[nzchar(x=var)]
  input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
  input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2)
  
  #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
  
  
  return(input_data)
}

# input_data <- generate_test_data()
# input_line <- generate_test_data_1()

# function to test/validate variant input
validate_input <- function(input_data) {
  
  input_data <- suppressWarnings(tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F))
  
  #chr = c(as.character(1:22),'x','X','y','Y')
  input_data <- input_data[ ( input_data$CHR.query %in% c(as.character(1:22),'x','X','y','Y')  &  grepl("^\\d", suppressWarnings(as.numeric(input_data$POS.query))) & input_data$REF.query %in% c("A","C", "T", "G") & input_data$ALT.query %in% c("A","C", "T", "G") ) , ]
  
}

# mb3 <- microbenchmark(out_new = (input_data$REF.query %in% c("A","C", "T", "G")),
#                       output_old = (grepl("[ATGC]", input_data$REF.query)  & nchar(input_data$REF.query)==1))
#  
# autoplot(mb3)

# input_line2 <- validate_input(input_line)
# 
# input_data2 <- validate_input(input_data)
# 


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

# input_data <- generate_test_data()
# input_line <- generate_test_data_1()
# input_data_2 <- generate_test_data_2()

# validate_input ----

# function to test/validate variant input
validate_input <- function(input_data) {
  
  input_data <- suppressWarnings(tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F))
  
  #chr = c(as.character(1:22),'x','X','y','Y')
  input_data <- input_data[ ( input_data$CHR.query %in% c(as.character(1:22),'x','X','y','Y')  &  grepl("^\\d", suppressWarnings(as.numeric(input_data$POS.query))) & input_data$REF.query %in% c("A","C", "T", "G") & input_data$ALT.query %in% c("A","C", "T", "G") ) , ]
  
}


# edit_paralog_colnames ----

paralog_colnames <- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  "ALT.query",
  "var.query",
  "ID.query",
  "Gene.query",
  "Codons.query",
  "Transcript.query",
  "Protein_dot.query",
  "Para_Z_score.query",
  "CHR.paralog",
  "POS.paralog",
  "REF.paralog",
  "ALT.paralog",
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
  "var.query",
  "Gene.query",
  "Positions.paralog")

paralog_extra_colnames <- c(
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


# lookup_vars ----

lookup_paralog <- function(input_data){
  
  paralog_out <- NULL
  #paralog_out <- read.csv(text = paste(paralog_colnames,collapse = ","))
  
  
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
      pg1 <- as.data.frame(stringr::str_split_fixed(tabix_paralog_extra, pattern = "\t", length(paralog_extra_colnames)), stringsAsFactors = F)
      
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
  colnames(paralog_out) <- paralog_extra_colnames
  
  #paralog_out <- paralog_out[,c(1:7,9:32)]
  #paralog_out <- dplyr::na_if(paralog_out, "NA")
  
  
  # if (nrow(paralog_out)!=0) {
  #   #paralog_out <- dplyr::na_if(paralog_out, "NA")
  #   return(paralog_out)
  # } else {
  #   paralog_out <- NA
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
  # if (nrow(paraloc_out)!=0) {
  #   #paralog_out <- dplyr::na_if(paralog_out, "NA")
  #   return(unique(paraloc_out))
  # } else {
  #   paraloc_out <- NA
  #   return(paraloc_out)
  # }
  return(unique(paraloc_out))
  
}

# predict_output ----



predict_output_tabix = function(input_data){


  #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F) 

    paralog <- lookup_paralog(input_data)
    paraloc <- lookup_paraloc(input_data)

  return(list("paralog" = paralog, "paraloc" = paraloc))

}  


# result <- predict_output_tabix(input_data = input_line)

  
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

# edit_download_colnames ----

# use tidyr::separate to split var.query and var.paralog columns for downloaded
# edit_download_cols = function(df) {
#   
#   #df <- tidyr::separate(df, var.query, into = c("CHR.query", "POS.query", "REF.query", "ALT.query") )
#   df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
#   
#   return(df)
# }

# edit_download_cols_paraloc = function(df) {
#   
#   df <- tidyr::separate(df, var, into = c("CHR.query", "POS.query", "REF.query") )
#   #df <- tidyr::separate(df, var.paralog, into = c("CHR.paralog", "POS.paralog", "REF.paralog", "ALT.paralog") )
#   colnames(df) <- c("CHR.query", "POS.query", "REF.query", "Gene.query", "Positions.paralog")
#   return(df)
# }

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

# function to add ensembl URL link
add_paraloc_gene_URLs <- function(pos) {
  
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
    
    
    # #Ensembl ENSG.query for paraloc
    # result_paraloc$ENSG.query <- ifelse(!is.na(result_paraloc$Gene.query),
    #                                     (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene.query)]), "' target='_blank'>", map[unlist(result_paraloc$Gene.query)], "</a>")),
    #                                     "-")
    
    #Ensembl Gene.query for paraloc
    result_paraloc$Gene.query<- ifelse(!is.na(result_paraloc$Gene.query), 
                                       (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene.query)]), "' target='_blank'>", result_paraloc$Gene.query, "</a>")),
                                       "-")
    
    
    
    for (i in c(1:nrow(result_paraloc))){ 
      # get all positions from list and apply add_URLs function to every position
      # then paste back as string
      result_paraloc$Positions.paralog[i] <- paste(unlist(lapply(unlist(result_paraloc$Positions.paralog[i]), function(line) add_paraloc_gene_URLs(line))), collapse = ", ")
      
    }
  } else { 
    result_paraloc <- NULL
  }
  


  return(result_paraloc)
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
    
    # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000213281
    #Ensembl ENSG.query
    # result$ENSG.query<- ifelse(!is.na(result$ENSG.query), 
    #                            (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.query), "' target='_blank'>", result$ENSG.query, "</a>")),
    #                            "-")
    # 
    #Ensembl ENSG.paralog
    # result$ENSG.paralog<- ifelse(!is.na(result$ENSG.paralog), 
    #                              (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.paralog), "' target='_blank'>", result$ENSG.paralog, "</a>")),
    #                              
    #                              "-")
    # 
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
             'Paralogue variant' = 'var.paralog',
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
  'Gene' = 'Gene.query',
  # 'ENSG' = 'ENSG.query',
  'Paralogous positions' = 'Positions.paralog')




# childrow_JS_callback ----

childrow_JS_callback <- c("
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
        '<tr>'+
            '<th>Query variant</th>'+
            '<td>'+ d[1] + '-' + d[2] + '-' + d[3] + '-' + d[4] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>ClinVar ID</th>'+
            '<td>'+ d[6] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>ClinVar Class</th>'+
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





####### draw proteins example ######

# "P51787 O43526 O43525 P56696 Q9NR82"

# data("five_rel_data")
# prot_data_example <- get_features("five_rel_data")
# # prot_data_example <- get_features("Q04206 Q01201 Q04864 P19838 Q00653")
# prot_data_example <- feature_to_dataframe(prot_data_example)


####################################

# example data
#input_data_2 <- generate_test_data_2(2) 
#result <- predict_output_tabix(validate_input(input_data_2))$paralog

# prepare draw protein dataframe ----

get_prot_data <- function(result) {
  
  # check number of input data vars
  if (length(unique(result$var.query))!=1) {
    result <- result[result$var.query == result$var.query[1],]
  }
  
  #paste(map_HGNC[unlist("KCNQ2")],collapse = " ")
  
  result$UniProt_ID.paralog <- map_HGNC[unlist(result$Gene.paralog)]
  result$UniProt_ID.query <- map_HGNC[unlist(result$Gene.query)]
  
  #query_proteins_rev <- paste(paste(rev(unique(result$UniProt_ID.paralog)),collapse = " "), map_HGNC[unlist(unique(result$Gene.query))],collapse = " ")
  
  query_proteins <- paste(map_HGNC[unlist(unique(result$Gene.query))], paste(unique(result$UniProt_ID.paralog),collapse = " "))
  query_HGVS <- paste(unique(result$ENST.query),"(", unique(result$Gene.query),"):",unique(result$cDNA.query), " (",unique(result$Protein.query),")", sep = "")
  
  prot_data <- feature_to_dataframe(get_features(query_proteins))
  
  
  # get pfam data
  pfam_data <- get_Pfam(query_proteins)
  # left join protein data
  pfam_data <- pfam_data %>% left_join(prot_data %>% select(accession, taxid, order) %>% unique.array(),
                                       by = "accession")
  pfam_data$begin <- as.numeric(pfam_data$begin)
  
  
  # jon Uniprot and Pfam data
  prot_data <- rbind(prot_data, pfam_data)
  
  
  # left join protein data
  prot_data <- prot_data %>% left_join(rbind(result %>% select("UniProt_ID" = UniProt_ID.query, "Gene" = Gene.query, "ENSG" = ENSG.query, "Protein_position" = Protein_position.query) %>% unique.array(),
                                             result %>% select("UniProt_ID" = UniProt_ID.paralog, "Gene" = Gene.paralog, "ENSG" = ENSG.paralog, "Protein_position" = Protein_position.paralog) %>% unique.array()),
                                       by = c("accession"= "UniProt_ID"))
    
  prot_data$Protein_position <- as.numeric(prot_data$Protein_position)
  prot_data$HGVS.query <- query_HGVS
  #prot_data$description <- 
  
  #Add Gene URL
  prot_data$Gene<- ifelse(
    prot_data$Gene!="NA",
    (paste0("<a href='", paste0("http://pfam.xfam.org/protein/",prot_data$accession,"/"), "' target='_blank'>", prot_data$Gene, "</a>")),
    prot_data$Gene)
  
  
  prot_data <- separate(data = prot_data,col = description, into = "description", extra = "drop", sep = ";")
  
  #list("prot_data" = prot_data, "query_HGVS" = query_HGVS)
  #return(list("prot_data" = prot_data, "query_HGVS" = query_HGVS))
  
  return(prot_data)
}
#####


# draw the prot graph ----
draw_prot_data_plotly <- function(input_data) {
  
  
  # get prot data
  prot_data <- get_prot_data( input_data)
  
  # uncomment and run for test data (1,2,3,4)
  #prot_data <- get_prot_data( predict_output_tabix(validate_input(generate_test_data_2(4)))$paralog)
  
  fig <- draw_plotly_graph(prot_data = prot_data)


  #fig_test <- subplot(draw_plotly_graph(prot_data = prot_data, showlegend = T), draw_plotly_graph(prot_data = prot_data_test, showlegend = F), nrows = 2, shareX = T) %>% layout(showlegend = T)


}

#####


# draw plotly graph ----

draw_plotly_graph <- function(prot_data, showlegend = T) {
  
  fig <- plot_ly(prot_data)
  
  # Chains
  if ("CHAIN" %in% prot_data$type ) {
    fig <- add_bars(fig, data = prot_data[prot_data$type == "CHAIN",],x = ~c((begin)-(end)), y = ~Gene, base = ~(end-Protein_position),
                    width = 0.02, orientation = 'h', showlegend = F , legendgroup = "Chains", # yaxis = ~Gene,  #xaxis = Gene
                    marker = list(color = toRGB("gray50")), name = ~Gene, hoverinfo = "name+text") # color = 'rgba(50, 171, 96, 0.6)'
  }
  
  
  # Folding
  if ("HELIX" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "HELIX",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.2, orientation = 'h', showlegend = F, name = ~type, legendgroup = "Folding",
                     opacity = 0.1 , hoverinfo = "name+text")
  }
  
  if ("STRAND" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "STRAND",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.2, orientation = 'h', showlegend = F, name = ~type,  legendgroup = "Folding",
                     opacity = 0.1 , hoverinfo = "name+text")
  }
  
  if ("TURN" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "TURN",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.2, orientation = 'h', showlegend = F, name = ~type,  legendgroup = "Folding",
                     opacity = 0.1 , hoverinfo = "name+text")
  }
  
  # Repeat
  if ("REPEAT" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "REPEAT",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.1, orientation = 'h', showlegend = F, name = ~description,  legendgroup = 'Repeat',
                     opacity = 0.2, marker = list(color = toRGB("gray50"),
                                                  line = list(color = toRGB("gray20"), width = 2)),
                     hoverinfo = "name+text")
  }
  
  # Region
  if ("REGION" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "REGION",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  # Domain
  if ("DOMAIN" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "DOMAIN",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  
  # Topo_Domain
  if ( ("TOPO_DOM" %in% prot_data$type ) | ("TRANSMEM" %in% prot_data$type)) {
    fig <- add_trace(fig, data = prot_data[(prot_data$type == "TOPO_DOM" | prot_data$type == "TRANSMEM") ,], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  
  
  # Motif
  if ("MOTIF" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "MOTIF",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'UniProt', 
                     hoverinfo = "name+text")
  }
  
  
  # PFam
  if ("PFAM" %in% prot_data$type ) {
    fig <- add_trace(fig, data = prot_data[prot_data$type == "PFAM",], x = ~c((begin)-(end)), y = ~Gene,type = 'bar', base = ~(end-Protein_position),
                     width = 0.4, orientation = 'h', showlegend = T, name = ~description,  legendgroup = 'PFam', 
                     hoverinfo = "name+text") 
  }
  
  # # Phosphorilation NOT WORKING
  # if ("MOD_RES" %in% prot_data$type) {
  #   fig <- add_trace(fig, data = phospho_site_info(prot_data), x = ~c((begin)-(end)), y = ~Gene,type = 'scatter', base = ~(end-Protein_position),
  #                    width = 0.04, orientation = 'h', showlegend = T, name = "Phosphorilation",  # legendgroup = 'PFam', 
  #                    hoverinfo = "name+text") 
  # }
  
  fig <- add_annotations(fig, data = prot_data,x = 0,y = 1, yref = "paper", yanchor = "bottom", showarrow = F, font = list(size = 14),
                         text = paste0(unique(prot_data$HGVS.query)) )  
  
  fig <- layout(fig, yaxis = list(title = "",autorange = T, showgrid = F, showline = F, showticklabels = T),  # , domain= c(0, 0.85) will set the total height of the plot 
                xaxis = list(title = "",autorange = T, showgrid = T, showline = F, showticklabels = F,zeroline = T, zerolinewidth = 3),
                barmode = 'overlay', showlegend = showlegend,legend = list(traceorder = "grouped", itemdoubleclick = "toggle", tracegroupgap = 20,
                                                   title = list(text = "Grouped Domain Annotations",font = list(size = 14)))) 
  fig <- config(fig, displayModeBar = T,
                modeBarButtonsToRemove = list("pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d","resetScale2d",
                                              "hoverClosestCartesian", "hoverCompareCartesian", "hoverClosestGl2d", "toggleSpikelines"),
                toImageButtonOptions = list(format = "png",
                                            width = "1800", height = "800"), displaylogo = F )
  
  
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
  
  if (nrow(pfam_eg)>=1) {
    print("Pfam download has worked.")
  }
  else {
    print(paste("An error has occured."))
  }
  
  #order=0
  pfam_eg[c("type", "description", "begin", "end", "length", "accession", "entryName")] <- NA
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
      
      #order <- order + 1
      #pfam_eg$order[i] <- order
      
    } else if (pfam_eg[i,1]=="A") {
      pfam_eg$type[i] <- "PFAM"
      pfam_eg$description[i] <- pfam_eg$X2[i]
      pfam_eg$begin[i] <- pfam_eg$X4[i]
      pfam_eg$end[i] <- pfam_eg$X5[i]
      pfam_eg$length[i] <- pfam_eg$X5[i] - as.numeric(pfam_eg$X4[i])
      pfam_eg$entryName[i] <- pfam_eg$X3[i]
      
      
      #pfam_eg$order[i] <- order 
      
    }
    if (is.na(pfam_eg$accession[i])){
      pfam_eg$accession[i] <- pfam_eg$accession[i-1]
    }
    
  }

  return(pfam_eg[,6:12])
}

######

# drawProteins functions ----
get_features <- function (proteins_acc) 
{
  proteins_acc_url <- gsub(" ", "%2C", proteins_acc)
  baseurl <- "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
  url <- paste0(baseurl, proteins_acc_url)
  prots_feat <- httr::GET(url, httr::accept_json())
  code <- httr::status_code(prots_feat)
  if (code == 200) {
    print("UniProt download has worked")
  }
  else {
    print(paste("An error has occured. Code:", code))
  }
  prots_feat_red <- httr::content(prots_feat)
  return(prots_feat_red)
}

feature_to_dataframe <- function (features_in_lists_of_six) 
{
  features_total_plot <- NULL
  for (i in 1:length(features_in_lists_of_six)) {
    features_temp <- extract_feat_acc(features_in_lists_of_six[[i]])
    features_temp$order <- i
    features_total_plot <- rbind(features_total_plot, features_temp)
  }
  return(features_total_plot)
}

extract_feat_acc <- function (features_list) 
{
  features <- NULL
  for (i in 1:length(features_list$features)) {
    if (is.null(features_list$features[[i]]$description) == 
        TRUE) {
      featuresTemp <- c(features_list$features[[i]]$type, 
                        "NONE", as.numeric(features_list$features[[i]]$begin), 
                        as.numeric(features_list$features[[i]]$end))
    }
    else {
      featuresTemp <- c(features_list$features[[i]]$type, 
                        as.character(features_list$features[[i]]$description), 
                        as.numeric(features_list$features[[i]]$begin), 
                        as.numeric(features_list$features[[i]]$end))
    }
    features <- rbind(features, featuresTemp)
  }
  features_dataframe <- as.data.frame(features, stringsAsFactors = FALSE)
  colnames(features_dataframe) <- c("type", "description", 
                                    "begin", "end")
  features_dataframe$begin <- as.numeric(features_dataframe$begin)
  features_dataframe$end <- as.numeric(features_dataframe$end)
  features_dataframe$length <- features_dataframe$end - features_dataframe$begin
  features_dataframe$accession <- rep(features_list$accession, 
                                      times = nrow(features_dataframe))
  features_dataframe$entryName <- rep(features_list$entryName, 
                                      times = nrow(features_dataframe))
  features_dataframe$taxid <- rep(features_list$taxid, times = nrow(features_dataframe))
  return(features_dataframe)
}

phospho_site_info <- function (features) 
{
  features <- features[features$type == "MOD_RES", ]
  phospho_list <- grep("Phospho", features$description)
  phospho_features <- features[phospho_list, ]
  return(phospho_features)
}



#######
