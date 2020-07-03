library(shiny)
library(DT)
library(shinythemes)
library(stringr)
library(tidyr)
# library(tidyverse)
# library(microbenchmark)


generate_test_data <- function() {
  
  # generate random input data for testing
  input1<-data.frame(chr="1",pos="115256528",ref="T",alt="C")
  input2<-data.frame(chr="1",pos="115256528",ref="T",alt="G")
  input3<-data.frame(chr="3",pos="38592567",ref="T",alt="A")
  input4<-data.frame(chr="21",pos="44592214",ref="C",alt="T")
  input5<-data.frame(chr="21",pos="47421902",ref="G",alt="A")
  input6<-data.frame(chr="X",pos="70443591",ref="G",alt="A")

  
  input <- rbind(input1,input2,input3,input4,input5,input6)
  
  #var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
  #var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
  var=var[nzchar(x=var)]
  input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
  input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
  input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2)
  
  return(input_data)
}

input_data <- generate_test_data()


#read gene symbol/ENSG and write to dict
# mart_export <- read.delim(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/mart_export.txt"), quote="", stringsAsFactors=F)
mart_export <- read.delim("data/mart_export.txt", quote="", stringsAsFactors=F)
map=setNames(mart_export$Gene.stable.ID, mart_export$HGNC.symbol)

HGNC_export <- read.delim("data/HGNC_all_genes.txt", quote="", stringsAsFactors=F)
#HGNC_export <- tidyr::separate(HGNC_export,HGNC_ID, into = c("HGNC", "ID"), remove = T)
map_HGNC = setNames(HGNC_export$HGNC_ID, HGNC_export$Approved_symbol)

# Load all data as Rds data
#raw_data = readRDS("./data/raw_data_paralog.Rds")
#Paraloc_data = readRDS("../../Paraloc_data_paraloc.Rds")


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

lookup_paralog <- function(input_data){
  
  paralog_out <- NULL
  #paralog_out <- read.csv(text = paste(paralog_colnames,collapse = ","))
  
  
  for (i in 1:nrow(input_data)) {
    
    #i=1
    query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
    #CMD_paralog<- paste0("tabix ", paralog_data, " ", query)
    #tabix_paralog <- system(command = paste0("tabix ", paralog_data, " ", query), intern = T,wait = T)

    #tabix_paralog <- system(command = paste0("tabix data/paralog_data_sorted.txt.gz ", query), intern = T,wait = T)
    tabix_paralog_extra <- system(command = paste0("tabix data/paralog_data_sorted.txt.gz ", query), intern = T,wait = T)
    
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
  }
  #paralog_out = rbind(paralog_out, pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),])
  
  #colnames(paralog_out) <- paralog_colnames
  colnames(paralog_out) <- paralog_extra_colnames
  
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
    tabix_paraloc <- system(command =  paste0("tabix data/paraloc_data_sorted.txt.gz " , query), intern = T, wait = T)
    
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
  paraloc_out$var.query <- stringr::str_replace_all(paraloc_out$var.query," ","-")
  paraloc_out$ENSG.query <- map[unlist(paraloc_out$Gene.query)]
  paraloc_out <- paraloc_out[,c(1:5,7,6)]
  
  return(unique(paraloc_out))
  
}




predict_output_tabix = function(input_data){


  input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F) 

  
  
  # call lookup function 
  # paraloc_out <- lookup_paraloc(input_data)
  
  return(list("output" = lookup_paralog(input_data), "paraloc_output" = lookup_paraloc(input_data)))
  
}  


#result <- predict_output_tabix(input_data = input_data)

  
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
    
    
    #Ensembl ENSG.query for paraloc
    result_paraloc$ENSG.query <- ifelse(!is.na(result_paraloc$Gene.query),
                                        (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",map[unlist(result_paraloc$Gene.query)]), "' target='_blank'>", map[unlist(result_paraloc$Gene.query)], "</a>")),
                                        "-")
    
    #Ensembl Gene.query for paraloc
    result_paraloc$Gene.query<- ifelse(!is.na(result_paraloc$Gene.query), 
                                       (paste0("<a href='", paste0("https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",map_HGNC[unlist(result_paraloc$Gene.query)]), "' target='_blank'>", result_paraloc$Gene.query, "</a>")),
                                       "-")
    
    
    
    for (i in c(1:nrow(result_paraloc))){ 
      # get all positions from list and apply add_URLs function to every position
      # then paste back as string
      result_paraloc$Positions.paralog[i] <- paste(unlist(lapply(unlist(result_paraloc$Positions.paralog[i]), function(line) add_URLs(line))), collapse = ", ")
      
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
    
    # https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000213281
    #Ensembl ENSG.query
    result$ENSG.query<- ifelse(!is.na(result$ENSG.query), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.query), "' target='_blank'>", result$ENSG.query, "</a>")),
                               "-")
    
    #Ensembl ENSG.paralog
    result$ENSG.paralog<- ifelse(!is.na(result$ENSG.paralog), 
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",result$ENSG.paralog), "' target='_blank'>", result$ENSG.paralog, "</a>")),
                                 
                                 "-")
    
    #HGNC Gene.query
    result$Gene.query<- ifelse(!is.na(result$Gene.query), 
                               (paste0("<a href='", paste0("https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",map_HGNC[unlist(result$Gene.query)]), "' target='_blank'>", result$Gene.query, "</a>")),
                               "-")
    
    #HGNC Gene.paralog
    result$Gene.paralog<- ifelse(!is.na(result$Gene.paralog), 
                                 (paste0("<a href='", paste0("https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/",map_HGNC[unlist(result$Gene.paralog)]), "' target='_blank'>", result$Gene.paralog, "</a>")),
                                 "-")
    
    
    return(result)
    
} 


paralog_DT_colnames <- c('Chr' = 'CHR.query',
             'Pos' = 'POS.query',
             'REF' = 'REF.query',
             'ALT' = 'ALT.query',
             'Query Variant' = 'var.query',
             'ClinVar ID' = 'ID.query',
             'Gene' = 'Gene.query',
             'ENSG' = 'ENSG.query',
             'ENST' = 'ENST.query',
             'cDNA' = 'cDNA.query',
             'Protein' = 'Protein.query',
             'cDNA position' = 'cDNA_position.query',
             'Protein position' = 'Protein_position.query',
             'AA' = 'AA.query',
             'Codons' = 'Codons.query',
             'Para_Z Score'='Para_Z_score.query',
             'Chr' = 'CHR.paralog',
             'Pos' = 'POS.paralog',
             'REF' = 'REF.paralog',
             'ALT' = 'ALT.paralog',
             'Paralogue variant' = 'var.paralog',
             'ClinVar ID' = 'ID.paralog',
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
  'ENSG' = 'ENSG.query',
  'Paralogous positions' = 'Positions.paralog')


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
    //        '<td>'+ d[9] + '(' + d[8] + '):' + d[10] + ' (' + d[11] + ')' +'</td>'+
    //        '<td>'+ d[15] +'</td>'+
    //        '<td>'+ d[16] +'</td>'+
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
            '<th>Gene</th>'+
            '<td>'+ d[7] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>HGVS</th>'+
            '<td>'+ d[9] + '(' + d[8] + '):' + d[10] + ' (' + d[11] + ')' +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Codons</th>'+
            '<td>'+ d[15] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Para_Z Score</th>'+
            '<td>'+ d[16] +'</td>'+
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



