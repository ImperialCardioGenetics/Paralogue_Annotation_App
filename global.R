library(shiny)
library(DT)
library(shinythemes)
library(tidyverse)
library(grid)
#library(arrow)

# max variants to keep from upload file
max_upload=50

# read clinvar
clinvar <- read_rds("data/clinvar_20220624.rds")

#read gene symbol/ENSG and write to dict
# HGNC <- read_feather("data/HGNC_complete_2022.feather")
HGNC <- read_rds("data/HGNC_complete_2022.rds")

# map_UniProt = setNames(HGNC$UniProt_ID, HGNC$Approved_symbol)
map_ENSG=setNames(HGNC$Ensembl_ID, HGNC$Approved_symbol)


# generate_test_data ----

# generate_test_data <- function() {
# 
#   # generate random input data for testing
#   input<-data.frame(chr="1",pos="115256528",ref="T")
#   input1<-data.frame(chr="1",pos="115256528",ref="T",alt="C")
#   input2<-data.frame(chr="1",pos="115256528DD",ref="T",alt="G")
#   input3<-data.frame(chr="3",pos="38592567",ref="TA",alt="B")
#   input4<-data.frame(chr="21",pos="44592214",ref="WC",alt="T")
#   input5<-data.frame(chr="21B",pos="47421902",ref="G",alt="A")
#   input6<-data.frame(chr="XB",pos="70443591",ref="G",alt="A")
# 
# 
#   input <- rbind(input1,input2,input3,input4,input5,input6)
# 
#   #var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
#   var = paste(input$chr,input$pos,input$ref,input$alt,sep = " ")
#   #var<-unlist(strsplit(input$var,split="\\, |\\,|\\n"))
#   var=var[nzchar(x=var)]
#   input_data<-data.frame(mutation=var, stringsAsFactors = FALSE)
#   input_data$mutation = stringr::str_replace_all(input_data$mutation,"[[:punct:][:space:]]","-")
#   input_data$mutation = stringr::str_replace_all(input_data$mutation,"^chr","")
#   input_data$paraloc = substr(input_data$mutation, 1, nchar(input_data$mutation)-2)
# 
#   #input_data <- tidyr::separate(input_data,mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F)
# 
# 
#   return(input_data)
# }

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

  } else if (var==6){
    #12-2690777-C-A  
    input<-data.frame(chr="12",pos="2690777",ref="C",alt="A")    

  } else if (var==7){
    #12-2690777-C-A  
    input0<-data.frame(chr="1",pos="115256513",ref="T",alt="G")
    input1<-data.frame(chr="12",pos="2690777",ref="C",alt="A")
    input2<-data.frame(chr="1",pos="115256528",ref="T",alt="G")
    input3<-data.frame(chr="1",pos="115256527",ref="T",alt="G")
    input4<-data.frame(chr="1",pos="24401872",ref="C",alt="A")
    input5<-data.frame(chr="1",pos="120572581",ref="A",alt="C")
    
    input <- rbind(input0,input1,input2,input3,input4,input5)
  
  } else if (var==8){
    #12-2690777-C-A  
    input<-data.frame(chr="1",pos="115256527",ref="T",alt="G")    
    
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
  
  if (!is.data.frame(input_data)) {
    input_data <- data.frame(mutation=unlist(strsplit(input_data,split="\\, |\\,|\\n|\\s|\\t")), stringsAsFactors = FALSE)
  } 
  
  # validate and filter any worng input vars
  input_data <- input_data %>%
    mutate(mutation = str_replace_all(mutation, "[[:punct:][:space:]]","-")) %>%
    mutate(mutation = str_replace_all(mutation, "^chr","")) %>%
    mutate(paraloc = str_sub(mutation, 1, nchar(mutation)-2)) %>%
    separate(mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F) %>%
    filter((CHR.query %in% c(as.character(1:22),'x','X','y','Y')) &
             (grepl("^\\d", suppressWarnings(as.numeric(POS.query)))) &
             (REF.query %in% c("A","C", "T", "G")) &
             (ALT.query %in% c("A","C", "T", "G")))
  
  
  return(input_data)
 
}



# check_upload_file ----

# function to check if uploaded variants file is valid
check_upload_file = function(inFile) {
  # very hacky way to read in vcf
  # read 1st line only to get number of cols to check if its txt or vcf format
  
  # inFile <- NA
  # inFile$datapath <- "data/test_data/test_upload1.txt"
  # inFile$datapath <- "data/test_data/test_upload2.txt"
  # inFile$datapath <- "data/test_data/test_upload3.txt"
  # inFile$datapath <- "data/test_data/chr19.txt"

  # inFile$datapath <- "data/test_data/chr1_head200.vcf"
  
  # check if upload file is gz
  if (endsWith(inFile$datapath, ".gz")) {
    input_1row = read.table(gzfile(inFile$datapath),nrows = 1,comment.char = "",colClasses = "character")
  } else {
    input_1row = read.table(inFile$datapath,nrows = 1,comment.char = "",colClasses = "character")
    #input_1row
  }
  # input_1row = read_delim(inFile$datapath,comment = "#",col_names = NA)
  # input_1row = read.table(file,nrows = 1,comment.char = "")
  
  # check if file is vcf
  if (str_starts(input_1row[1,1],"##fileformat=VCF") | str_starts(input_1row[1,1],"#CHROM")) {
    
    if (ncol(input_1row)==1 | ncol(input_1row)>=5) {
      #input_file <- read_delim(inFile$datapath, delim = "\t",col_select = (if_else ((ncol(input_1row)==1 | ncol(input_1row)>=5),as_vector(c(1,2,4,5)),as_vector(c(1,2,3,4)))),comment = "#",col_names = NA,show_col_types = F,col_types = cols(.default = "c")) %>%
      #input_file <- read_delim(inFile$datapath, delim = "\t",col_select = (if_else ((ncol(input_1row)==1 | ncol(input_1row)>=5),select(c(1,2,4,5)),select(c(1,2,3,4)))),comment = "#",col_names = NA,show_col_types = F,col_types = cols(.default = "c")) %>%
      input_file <- read_delim(inFile$datapath, delim = "\t",col_select = c(1,2,4,5),comment = "#",col_names = NA,show_col_types = F,col_types = cols(.default = "c")) %>%
        slice_head(n=max_upload) %>%  # Keep only first 50 vars
        #do(if(ncol(.) <= 4) . else select(.,c(1,2,4,5))) %>% 
        filter(nchar(.[[3]])==nchar(.[[4]])) %>% 
        unite(sep = ":",remove = T,col = "mutation") %>% 
        filter((nchar(mutation)<10) | (str_ends(mutation,"[AGTC]$")) | (str_detect(mutation,"\\*")))
      #input_file <- read_delim("data/test_data/chr1_head200.vcf", delim = "\t",col_select = c(1,2,4,5),comment = "#",col_names = NA) %>% unite(sep = ":",remove = T,col = "col1")
    } else if (ncol(input_1row)==4) {
      input_file <- read_delim(inFile$datapath, delim = "\t",col_select = c(1,2,3,4),comment = "#",col_names = NA,show_col_types = F,col_types = cols(.default = "c")) %>%
        slice_head(n=max_upload) %>%   # Keep only first 50 vars
        #do(if(ncol(.) <= 4) . else select(.,c(1,2,4,5))) %>% 
        filter(nchar(.[[3]])==nchar(.[[4]])) %>% 
        unite(sep = ":",remove = T,col = "mutation") %>% 
        filter((nchar(mutation)<10) | (str_ends(mutation,"[AGTC]$")) | (str_detect(mutation,"\\*")))
    }
      #input_file <- read_delim("data/test_data/chr1_head200.vcf", delim = "\t",col_select = c(1,2,4,5),comment = "#",col_names = NA) %>% unite(sep = ":",remove = T,col = "col1")
    # check if file is flat txt
  } else if (ncol(input_1row)>=1 & ncol(input_1row)<=4 &  grepl("^chr|^[1-9]|^[XY]", input_1row[1])) {
    
    #input_file <- read_delim(inFile$datapath, delim = "[:punct:]", comment = "##",col_names = NA, show_col_types = F,col_types = cols(.default = "c")) #%>% 
    input_file <- read.table(text = gsub("[[:punct:][:space:]]","\t",read_lines(inFile$datapath) %>% str_subset(pattern = "\\*",negate = T))) %>% 
      slice_head(n=max_upload) %>%  # Keep only first 50 vars
      filter(nchar(.[[3]])==nchar(.[[4]])) %>%
      unite(sep = ":",remove = T,col = "mutation") %>% 
      filter((nchar(mutation)<10) | (str_ends(mutation,"[AGTC]$")))
    
  } else {
    # write NA table
    input_file<- tibble(mutation=NA)
  }
  
 

  # print("file uploaded")
  # print(paste0("n vars ",nrow(input_file)))
  return(input_file)
}



# edit_paralog_colnames ----

paraloc_colnames<- c(
  "CHR.query",
  "POS.query",
  "REF.query",
  #"var.query",
  "Gene.query",
  "Positions.paraloc")

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

      #tabix_paralog <- system(command = paste0("tabix data/paralog_data_sorted.txt.gz ", query), intern = T,wait = T)
      tabix_paralog_extra <- suppressMessages(suppressWarnings(system(command = paste0("tabix data/paralog_data.txt.gz ", query), intern = T,wait = T)))

      pg1 <- as.data.frame(stringr::str_split_fixed(tabix_paralog_extra, pattern = "\t", length(paralog_colnames)), stringsAsFactors = F)
      
      #colnames(pg1) <- paralog_colnames
      
      
      if (is.null(paralog_out)){
        paralog_out = pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),]
      } else {
        paralog_out = rbind(paralog_out, pg1[(pg1$V3==input_data$REF.query[i] & pg1$V4 == input_data$ALT.query[i]),])
      }
    #}
  }

  colnames(paralog_out) <- paralog_colnames
  
  # sort df to avoid sorting later
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
    
    tabix_paraloc <-suppressMessages(suppressWarnings(system(command =  paste0("tabix data/paraloc_chr/paraloc_data_chr", input_data$CHR.query[i] ,".txt.gz " , query), intern = T, wait = T)))
    
    pc1<- as.data.frame(stringr::str_split_fixed(tabix_paraloc, pattern = "\t", n = length(paraloc_colnames)), stringsAsFactors = F)

    if (is.null(paraloc_out)){
      paraloc_out = pc1[(pc1$V3==input_data$REF.query[i]),]
    } else {
      paraloc_out = rbind(paraloc_out, pc1[(pc1$V3==input_data$REF.query[i]),])
    }
  }
  
  colnames(paraloc_out) <- paraloc_colnames

  return(unique(paraloc_out))
  
}


lookup_homolog <- function(input_data){
  
  homolog_out <- NULL

  for (i in 1:nrow(input_data)) {
    
    #i=1
    query <- paste0(input_data$CHR.query[i], ":", input_data$POS.query[i], "-", input_data$POS.query[i])
    tabix_homolog <- suppressMessages(suppressWarnings(system(command = paste0("tabix data/homolog_data.txt.gz ", query), intern = T,wait = T)))

    hg1 <- as.data.frame(stringr::str_split_fixed(tabix_homolog, pattern = "\t", length(homolog_colnames)), stringsAsFactors = F)

    if (is.null(homolog_out)){
      homolog_out = hg1[(hg1$V3==input_data$REF.query[i] & hg1$V4 == input_data$ALT.query[i]),]
    } else {
      homolog_out = rbind(homolog_out, hg1[(hg1$V3==input_data$REF.query[i] & hg1$V4 == input_data$ALT.query[i]),])
    }
  }
  
  colnames(homolog_out) <- homolog_colnames
  
  # sort df to avoid sorting later
  homolog_out <- homolog_out[order(factor(homolog_out$CHR.query , levels = c(1:22,"X","Y")), 
                                   as.numeric(homolog_out$POS.query),
                                   factor(homolog_out$CHR.homolog , levels = c(1:22,"X","Y")), 
                                   as.numeric(homolog_out$POS.homolog)), ]

  return(homolog_out)
}


# lookup paralog new
lookup_paralog_new <- function(input_data) {
  
  paralog <- input_data %>% 
    left_join(clinvar,by = c("CHR.query" = "CHROM.c" ,"POS.query" = "POS.c", "REF.query" = "REF.c"),keep=T) %>% 
    rename_with(~ gsub(.,pattern = "\\.c",replacement = "\\.query.c")) %>% 
    rowwise() %>% 
    mutate(tbx_out = list(suppressMessages(suppressWarnings(system(paste0("tabix data/paraloc_chr/paraloc_data_chr",CHR.query,".txt.gz ",CHR.query,":",POS.query,"-",POS.query), intern = T))))) %>% 
    separate(tbx_out, into = c("CHR.query.p","POS.query.p","REF.query.p","Gene.query.p","Positions.paraloc.p"),sep = "\t",fill = "right")  %>%
    filter(REF.query==REF.query.p) %>% 
    separate_rows(Positions.paraloc.p,sep = ", ") %>%
    separate(Positions.paraloc.p,into = c("Gene.paraloc", "CHR.paraloc", "AA_pos.paraloc","AA.paraloc"),sep = " ") %>%
    separate(Gene.paraloc, into = c("Gene.paraloc","Gene.ext"),"-",fill = "right") %>% 
    select(-Gene.ext) %>% 
    separate(AA_pos.paraloc, into = c("AA_start.paraloc","AA_end.paraloc"),remove = F) %>%
    mutate(POS.paraloc = map2(AA_start.paraloc, AA_end.paraloc, seq)) %>%
    unnest(POS.paraloc) %>% 
    mutate(across(everything(), as.character)) %>% 
    left_join(clinvar,by = c("CHR.paraloc" = "CHROM.c" ,"POS.paraloc" = "POS.c"),keep=T) %>% 
    filter(!is.na(CHROM.query.c), !is.na(CHROM.c),ALT.query==ALT.query.c) %>% 
    select("CHR.query"               = CHROM.query.c,
           "POS.query"               = POS.query.c,
           "REF.query"               = REF.query.c,
           "ALT.query"               = ALT.query.c,
           "var_id.query"            = var_id.query.c,
           "ID.query"                = ID.query.c,
           "CLNSIG.query"            = CLNSIG.query.c,
           "SYMBOL.query"            = SYMBOL.query.c,
           "Gene.query"              = Gene.query.c,
           "Feature.query"           = Feature.query.c,
           "cDNA.query"              = cDNA.query.c,
           "Amino_acids.query"       = Amino_acids.query.c,
           "cDNA_position.query"     = cDNA_position.query.c,
           "AA_position.query"       = Protein_position.query.c,
           "Codons.query"            = Codons.query.c,
           "Para_z_score.query"      = Para_z_score.query.c,
           "CHR.paralog"             = CHROM.c,
           "POS.paralog"             = POS.c,
           "REF.paralog"             = REF.c,
           "ALT.paralog"             = ALT.c,
           "var_id.paralog"          = var_id.c,
           "ID.paralog"              = ID.c,
           "CLNSIG.paralog"          = CLNSIG.c,
           "SYMBOL.paralog"          = SYMBOL.c,
           "Gene.paralog"            = Gene.c,
           "Feature.paralog"         = Feature.c,
           "cDNA.paralog"            = cDNA.c,
           "Amino_acids.paralog"     = Amino_acids.c,
           "cDNA_position.paralog"   = cDNA_position.c,
           "AA_position.paralog"     = Protein_position.c,
           "Codons.paralog"          = Codons.c,
           "Para_z_score.paralog"    = Para_z_score.c) %>% 
    distinct() %>%
    arrange(factor(CHR.query , levels = c(1:22,"X","Y")),
            as.numeric(POS.query),
            factor(CHR.paralog , levels = c(1:22,"X","Y")),
            as.numeric(POS.query))
  
  return(paralog)
  
}


# lookup paraloc new
lookup_paraloc_new <- function(input_data) {
  
  paraloc <- input_data %>% 
    left_join(clinvar,by = c("CHR.query" = "CHROM.c" ,"POS.query" = "POS.c", "REF.query" = "REF.c"),keep=T) %>% 
    rename_with(~ gsub(.,pattern = "\\.c",replacement = "\\.query.c")) %>% 
    rowwise() %>% 
    mutate(tbx_out = list(suppressMessages(suppressWarnings(system(paste0("tabix data/paraloc_chr/paraloc_data_chr",CHR.query,".txt.gz ",CHR.query,":",POS.query,"-",POS.query), intern = T))))) %>% 
    separate(tbx_out, into = c("CHR.query.p","POS.query.p","REF.query.p","Gene.query.p","Positions.paraloc.p"),sep = "\t",fill = "right")  %>%
    filter(REF.query==REF.query.p) %>% 
    separate_rows(Positions.paraloc.p,sep = ", ") %>%
    separate(Positions.paraloc.p,into = c("Gene.paraloc", "CHR.paraloc", "AA_pos.paraloc","AA.paraloc"),sep = " ") %>%
    separate(Gene.paraloc, into = c("Gene.paraloc","Gene.ext"),"-",fill = "right") %>% 
    select(-Gene.ext) %>% 
    separate(AA_pos.paraloc, into = c("AA_start.paraloc","AA_end.paraloc"),remove = F) %>%
    mutate(POS.paraloc = map2(AA_start.paraloc, AA_end.paraloc, seq)) %>%
    unnest(POS.paraloc) %>% 
    mutate(across(everything(), as.character)) %>% 
    left_join(clinvar,by = c("CHR.paraloc" = "CHROM.c" ,"POS.paraloc" = "POS.c"),keep=T) %>% 
    select("CHR.query"               = CHR.query.p,
             "POS.query"               = POS.query.p,
             "REF.query"               = REF.query.p,
             "var_id.query"            = paraloc,
             "Gene.query"              = Gene.query.p,
             "AA.query"                = AA.paraloc,
             "Gene.paraloc"            = Gene.paraloc,
             "CHR.paraloc"             = CHR.paraloc,
             "AA_pos.paraloc"          = AA_pos.paraloc) %>% 
    distinct() %>% 
    arrange(factor(CHR.query , levels = c(1:22,"X","Y")),
            as.numeric(POS.query),
            factor(CHR.paraloc , levels = c(1:22,"X","Y"))) #%>% as_tibble()
 return(paraloc)

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
add_paraloc_URL_new = function(paraloc) {
  #if (nrow(paraloc)!=0){
    #Ensembl Gene.query for paraloc
  
  paraloc$Gene_ENSG.query   <- map_ENSG[unlist(paraloc$Gene.query)] 
  paraloc$Gene_ENSG.paraloc <- map_ENSG[unlist(paraloc$Gene.paraloc)] 

    paraloc$Gene.query <- ifelse(paraloc$Gene.query!="NA",
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",paraloc$Gene_ENSG.query), "' target='_blank'>", paraloc$Gene.query, "</a>")),
                                 "-")
    #Ensembl gene for paraloc
    paraloc$Gene.paraloc<- ifelse(paraloc$Gene.paraloc!="NA",
                                  (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",paraloc$Gene_ENSG.paraloc), "' target='_blank'>", paraloc$Gene.paraloc, "</a>")),
                                  "-")
    paraloc <- select(paraloc,-c("Gene_ENSG.query","Gene_ENSG.paraloc"))

}

# function to add ensembl URLs to paralg positions in genes
add_paraloc_URL_old = function(result_paraloc) {

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
add_paralog_URL = function(paralog) {
  

    #ClinVarID paralog URL
    paralog$ID.paralog<- ifelse(#!is.na(paralog$ID.paralog),
      paralog$ID.paralog!="NA",
      (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",paralog$ID.paralog,"/"), "' target='_blank'>", paralog$ID.paralog, "</a>")),
      "-")
    
    #ClinVarID query URL
    paralog$ID.query<- ifelse(#!is.na(paralog$ID.query),
      paralog$ID.query!="NA",
      (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",paralog$ID.query,"/"), "' target='_blank'>", paralog$ID.query, "</a>")),
      "-")
    #print(paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(paralog$Gene.query)],";g1=",map[unlist(paralog$SYMBOL.paralog)]))
    #Ensembl alignment URL
    # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
    paralog$Ensembl_alignment_link<- ifelse(!is.na(paralog$Gene.paralog), 
                                           (paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",paralog$Gene.query,";g1=",paralog$Gene.paralog), "' class='btn btn-default btn-xs btn-block active' target='_blank'>alignment</a>")) , 
                                           "-") 
    
    # https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000213281;t=ENST00000369535
    #Ensembl ENST.query
    paralog$Feature.query<- ifelse(!is.na(paralog$Feature.query), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",paralog$Gene.query,";t=",paralog$Feature.query), "' target='_blank'>", paralog$Feature.query, "</a>")),
                               "-")
    
    #Ensembl ENST.paralog
    paralog$Feature.paralog<- ifelse(!is.na(paralog$Feature.paralog), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",paralog$Gene.paralog,";t=",paralog$Feature.paralog), "' target='_blank'>", paralog$Feature.paralog, "</a>")),
                               "-")
    
    #HGNC Gene.query
    paralog$SYMBOL.query<- ifelse(!is.na(paralog$SYMBOL.query), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",paralog$Gene.query), "' target='_blank'>", paralog$SYMBOL.query, "</a>")),
                               "-")
    
    #HGNC Gene.paralog
    paralog$SYMBOL.paralog<- ifelse(!is.na(paralog$SYMBOL.paralog), 
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",paralog$Gene.paralog), "' target='_blank'>", paralog$SYMBOL.paralog, "</a>")),
                                 "-")
    
    paralog <- cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', paralog )
    
    return(paralog)
    
} 

# function to add ensembl URLs to homolog positions in genes
add_homolog_URL = function(result) {
  
  
  result$Pfam_pos.query <- strsplit(result$Pfam_pos.query,split = "=")[[1]][2]
  # https://pfam.xfam.org/family/
  
  #ClinVarID paralog URL
  result$Pfam_domain.query<- ifelse(#!is.na(result$ID.paralog),
    result$Pfam_domain.query!="NA",
    (paste0("<a href='", paste0("https://pfam.xfam.org/family/",result$Pfam_domain.query,"/"), "' target='_blank'>", result$Pfam_domain.query, "</a>")),
    "-")
  
  
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
                         'Query variant' = 'var_id.query',
                         'ClinVar.query' = 'ID.query',
                         'ClinVar Class.query' = 'CLNSIG.query',
                         'Gene.query' = 'SYMBOL.query',
                         'ENSG.query' = 'Gene.query',
                         'ENST.query' = 'Feature.query',
                         'cDNA.query' = 'cDNA.query',
                         'Residue.query' = 'Amino_acids.query',
                         'cDNA position.query' = 'cDNA_position.query',
                         'Residue_position.query' = 'AA_position.query',
                         'Codons.query' = 'Codons.query',
                         'Para_Z Score.query'='Para_z_score.query',
                         'Chr' = 'CHR.paralog',
                         'Pos' = 'POS.paralog',
                         'REF' = 'REF.paralog',
                         'ALT' = 'ALT.paralog',
                         'Paralogous variant' = 'var_id.paralog',
                         'ClinVar ID' = 'ID.paralog',
                         'ClinVar Class' = 'CLNSIG.paralog',
                         'Gene' = 'SYMBOL.paralog',
                         'ENSG' = 'Gene.paralog',
                         'ENST' = 'Feature.paralog',
                         'cDNA' = 'cDNA.paralog',
                         'Residue' = 'Amino_acids.paralog',
                         'cDNA position' = 'cDNA_position.paralog',
                         'Residue position' = 'AA_position.paralog',
                         'Codons' = 'Codons.paralog',
                         'Para_Z Score' = 'Para_z_score.paralog',
                         'Ensembl alignment' = 'Ensembl_alignment_link')



paraloc_DT_colnames <-  c('POS' = 'CHR.query', 
                          'CHR' = 'POS.query', 
                          'REF' = 'REF.query', 
                          'Query variant' = 'var_id.query',
                          'Query gene' = 'Gene.query',
                          'Query residue' = 'AA.query',
                          'Paralogous gene' = 'Gene.paraloc',
                          'Chromosome' = 'CHR.paraloc',
                          'AA positions' = 'AA_pos.paraloc'
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
                         'Residue.query' = 'Protein.query',
                         'cDNA position.query' = 'cDNA_position.query',
                         'Residue_position.query' = 'Protein_position.query',
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
                         'Residue' = 'Protein.homolog',
                         'cDNA position' = 'cDNA_position.homolog',
                         'Residue position' = 'Protein_position.homolog',
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
    //        '<td>'+ d[15] +'</td>'+
    //        '<td>'+ d[16] +'</td>'+
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


# DT options ----
paralog_DT_options_list = list(
  rowGroup = list(dataSrc = c(5)),
  #dom = 'lfrti',
  dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" +
         "<"row"<"col-sm-12"tr>>" +
         "<"row"<"col-sm-6"i>>" +
         "<"row"<"col-sm-6 btn-md"B>>"',
  buttons = list(list(extend = 'excel',
                      text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                      filename = paste0("paralogous_annotations_",Sys.Date())),
                 list(extend = 'csv',
                      fieldBoundary = '',
                      text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                      fieldSeparator = '\t',
                      filename = paste0("paralogous_annotations_",Sys.Date()),
                      extension = '.txt')),
  paging = T,scrollX = TRUE,
  columnDefs = list(
    list(visible = FALSE, targets = c(1:4,6:20,25,29,30)),
    list(width = "100px",targets = 5),
    list(orderable = FALSE, className = 'details-control', targets = 0)))


paraloc_DT_options_list = list(
  rowGroup = list(dataSrc = c(3)),
  #dom = 'lfrti',
  dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" + 
         "<"row"<"col-sm-12"tr>>" + 
         "<"row"<"col-sm-6"i>>" +
         "<"row"<"col-sm-6 btn-md"B>>"',
  buttons = list(list(extend = 'excel',
                      text = ' <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                      filename = paste0("paralogous_positions_",Sys.Date())),
                 list(extend = 'csv',
                      fieldBoundary = '',
                      text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                      fieldSeparator = '\t',
                      filename = paste0("paralogous_positions_",Sys.Date()),
                      extension = '.txt')),
  paging = T,scrollX = FALSE,
  columnDefs = list(
    list(visible = FALSE, targets = c(0:2)),#,9,10 )),# delete the 9,10 to add API calls,9,10 ENST, prot_pos)),
    list(width = "140px",targets = 3),
    list(className = 'dt-center', targets = c(3))))


homolog_DT_options_list = list(
  rowGroup = list(dataSrc = c(5)),
  #dom = 'lfrti',
  dom = '"<"row"<"col-sm-6"l><"col-sm-6"f>>" +
         "<"row"<"col-sm-12"tr>>" +
         "<"row"<"col-sm-6"i>>" +
         "<"row"<"col-sm-6 btn-md"B>>"',
  buttons = list(list(extend = 'excel',
                      text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.xslx)',
                      filename = paste0("homologous_pfam_annotations_",Sys.Date())),
                 list(extend = 'csv',
                      fieldBoundary = '',
                      text = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-download"></i>  Download (.txt)',
                      fieldSeparator = '\t',
                      filename = paste0("homologous_pfam_annotations_",Sys.Date()),
                      extension = '.txt')),
  paging = T,scrollX = TRUE,
  columnDefs = list(
    list(visible = FALSE, targets = c(1:4,6:22,27, 31:33)),
    list(width = "100px",targets = 5),
    list(orderable = FALSE, className = 'details-control', targets = 0)))

#######
