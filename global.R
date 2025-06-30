suppressWarnings(suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(tidyverse)
  library(DT)
}))

# max variants to keep from upload file
max_upload=50

# read clinvar
# load("data/clinvar_20241103_clean.RData")

clinvar_version="20250623"

load(paste0("data/clinvar_",clinvar_version,"_clean.RData"))
load("data/homolog_data_raw.RData")


# generate_test_data ----


generate_test_data <- function(var=NULL) {
  
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
    
  } else if (var==9){
    #12-2690777-C-A  
    input<-data.frame(chr="1",pos="237804283",ref="G",alt="A")  
    
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
    separate(mutation, into = c("CHR.query", "POS.query", "REF.query", "ALT.query"), remove = F,fill = "right") %>%
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





# lookup_vars ----


# lookup paralog new
lookup_paralog_new <- function(input_data) {
  
  #input_data <- validate_input(generate_test_data(9))
  
  paralog <- input_data %>% 
    left_join(clinvar,by = c("CHR.query" = "CHROM.c" ,"POS.query" = "POS.c", "REF.query" = "REF.c"),keep=T,relationship = "many-to-many") %>% 
    rename_with(~ gsub(.,pattern = "\\.c",replacement = "\\.query.c")) %>% 
    rowwise() %>% 
    mutate(tbx_out = list(suppressMessages(suppressWarnings(system(paste0("tabix data/paraloc_chr/paraloc_data_chr",CHR.query,".txt.gz ",CHR.query,":",POS.query,"-",POS.query), intern = T))))) %>% 
    separate(tbx_out, into = c("CHR.query.p","POS.query.p","REF.query.p","Gene.query.p","Positions.paraloc.p"),sep = "\t",fill = "right")  %>%
    filter(REF.query==REF.query.p) %>% 
    separate_rows(Positions.paraloc.p,sep = ", ") %>%
    separate(Positions.paraloc.p,into = c("Gene.paraloc", "CHR.paraloc", "AA_pos.paraloc","AA.paraloc"),sep = " ",fill = "right") %>%
    separate(Gene.paraloc, into = c("Gene.paraloc","Gene.ext"),"-",fill = "right") %>% 
    select(-Gene.ext) %>% 
    separate(AA_pos.paraloc, into = c("AA_start.paraloc","AA_end.paraloc"),remove = F,fill = "right") %>%
    mutate(POS.paraloc = map2(AA_start.paraloc, AA_end.paraloc, seq)) %>%
    unnest(POS.paraloc) %>% 
    mutate(across(everything(), as.character)) %>% 
    left_join(clinvar,by = c("CHR.paraloc" = "CHROM.c" ,"POS.paraloc" = "POS.c"),keep=T,relationship = "many-to-many") %>% 
    filter(!is.na(CHROM.query.c), !is.na(CHROM.c),ALT.query==ALT.query.c) %>% 
    separate(HGVSc.query.c, into = c("Feature.query.c","cdot.query.c"),sep = ":",remove = F,fill = "right") %>% 
    separate(HGVSc.c, into = c("Feature.c","cdot.c"),sep = ":",remove = F,fill = "right") %>% 
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
           "cdot.query"              = cdot.query.c,           
           "HGVSc.query"             = HGVSc.query.c,
           "HGVSp.query"             = HGVSp.query.c,
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
           "cdot.paralog"            = cdot.c,
           "HGVSc.paralog"           = HGVSc.c,
           "HGVSp.paralog"           = HGVSp.c,
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
    left_join(clinvar,by = c("CHR.query" = "CHROM.c" ,"POS.query" = "POS.c", "REF.query" = "REF.c"),keep=T,relationship = "many-to-many") %>% 
    rename_with(~ gsub(.,pattern = "\\.c",replacement = "\\.query.c")) %>% 
    rowwise() %>% 
    mutate(tbx_out = list(suppressMessages(suppressWarnings(system(paste0("tabix data/paraloc_chr/paraloc_data_chr",CHR.query,".txt.gz ",CHR.query,":",POS.query,"-",POS.query), intern = T))))) %>% 
    separate(tbx_out, into = c("CHR.query.p","POS.query.p","REF.query.p","Gene.query.p","Positions.paraloc.p"),sep = "\t",fill = "right")  %>%
    filter(REF.query==REF.query.p) %>% 
    separate_rows(Positions.paraloc.p,sep = ", ") %>%
    separate(Positions.paraloc.p,into = c("Gene.paraloc", "CHR.paraloc", "AA_pos.paraloc","AA.paraloc"),sep = " ",fill = "right") %>%
    separate(Gene.paraloc, into = c("Gene.paraloc","Gene.ext"),sep = "-",fill = "right") %>% 
    select(-Gene.ext) %>% 
    separate(AA_pos.paraloc, into = c("AA_start.paraloc","AA_end.paraloc"),remove = F,fill = "right") %>%
    mutate(POS.paraloc = map2(AA_start.paraloc, AA_end.paraloc, seq)) %>%
    unnest(POS.paraloc) %>% 
    mutate(across(everything(), as.character)) %>% 
    left_join(clinvar,by = c("CHR.paraloc" = "CHROM.c" ,"POS.paraloc" = "POS.c"),keep=T,relationship = "many-to-many") %>% 
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



# lookup homolog new
lookup_homolog_new <- function(input_data) {
  
  homolog <- input_data %>% 
    left_join(clinvar,by = c("CHR.query" = "CHROM.c" ,"POS.query" = "POS.c", "REF.query" = "REF.c"),keep=T,relationship = "many-to-many") %>% 
    rename_with(~ gsub(.,pattern = "\\.c",replacement = "\\.query.c")) %>%
    left_join(homolog_data_raw,by = c("CHROM.query.c" = "CHROM.query" ,"POS.query.c" = "POS.query", "REF.query.c" = "REF.query","ALT.query.c" = "ALT.query"),relationship = "many-to-many") %>% 
    left_join(clinvar,by = c("CHROM.homolog" = "CHROM.c" ,"POS.homolog" = "POS.c", "REF.homolog" = "REF.c","ALT.homolog" = "ALT.c"),keep=T,relationship = "many-to-many") %>% 
    filter(!is.na(var_id.query.c),!is.na(var_id.c),ALT.query==ALT.query.c) %>% 
    separate(HGVSc.query.c, into = c("Feature.query.c","cdot.query.c"),sep = ":",remove = F,fill = "right") %>% 
    separate(HGVSc.c, into = c("Feature.c","cdot.c"),sep = ":",remove = F,fill = "right") %>% 
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
           "cdot.query"              = cdot.query.c,           
           "HGVSc.query"             = HGVSc.query.c,
           "HGVSp.query"             = HGVSp.query.c,
           "Amino_acids.query"       = Amino_acids.query.c,
           "cDNA_position.query"     = cDNA_position.query.c,
           "AA_position.query"       = Protein_position.query.c,
           "Codons.query"            = Codons.query.c,
           "Pfam_domain.query"       = Pfam_domain.query,
           "Pfam_pos.query"          = Pfam_pos.query,
           "CHR.homolog"             = CHROM.c,
           "POS.homolog"             = POS.c,
           "REF.homolog"             = REF.c,
           "ALT.homolog"             = ALT.c,
           "var_id.homolog"          = var_id.c,
           "ID.homolog"              = ID.c,
           "CLNSIG.homolog"          = CLNSIG.c,
           "SYMBOL.homolog"          = SYMBOL.c,
           "Gene.homolog"            = Gene.c,
           "Feature.homolog"         = Feature.c,
           "cdot.homolog"            = cdot.c,
           "HGVSc.homolog"           = HGVSc.c,
           "HGVSp.homolog"           = HGVSp.c,
           "Amino_acids.homolog"     = Amino_acids.c,
           "cDNA_position.homolog"   = cDNA_position.c,
           "AA_position.homolog"     = Protein_position.c,
           "Codons.homolog"          = Codons.c) %>% 
    distinct() %>%
    arrange(factor(CHR.query , levels = c(1:22,"X","Y")),
            as.numeric(POS.query),
            factor(CHR.homolog , levels = c(1:22,"X","Y")),
            as.numeric(POS.query))
  
  return(homolog)
  
}



# test_lookup_functions ----

# test functions with random var
# result_paralog <- lookup_paralog_new(validate_input(generate_test_data(9)))
# result_paraloc <- lookup_paraloc_new(validate_input(generate_test_data(1)))
# result_homolog <- lookup_homolog_new(validate_input(generate_test_data(1)))



# add_URLs ----

# function to add ensembl URLs to positions NEW
add_URL = function(data) {
  
  if (all(c("Pfam_domain.query","Pfam_pos.query") %in% names(data))) {
    
    # data$Pfam_pos.query <- strsplit(data$Pfam_pos.query,split = "=")[[1]][2]
    # https://pfam.xfam.org/family/
    
    #ClinVarID homolog URL
    data$Pfam_domain.query<- ifelse(data$Pfam_domain.query!="NA",
                                    (paste0("<a href='", paste0("https://pfam.xfam.org/family/",data$Pfam_domain.query,"/"), "' target='_blank'>", data$Pfam_domain.query, "</a>")),
                                    "-")
    
    #ClinVarID ID.query
    data$ID.query <- ifelse(data$ID.query!="NA",
                            (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",data$ID.query,"/"), "' target='_blank'>", data$ID.query, "</a>")),
                            "-")
    
    #ClinVarID ID.homolog
    data$ID.homolog<- ifelse(data$ID.homolog!="NA",
                             (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",data$ID.homolog,"/"), "' target='_blank'>", data$ID.homolog, "</a>")),
                             "-")
    
    #Ensembl SYMBOL.query
    data$SYMBOL.query <- ifelse(!is.na(data$SYMBOL.query), 
                                (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",data$Gene.query), "' target='_blank'>", data$SYMBOL.query, "</a>")),
                                "-")
    
    #Ensembl SYMBOL.homolog
    data$SYMBOL.homolog<- ifelse(!is.na(data$SYMBOL.homolog), 
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",data$Gene.homolog), "' target='_blank'>", data$SYMBOL.homolog, "</a>")),
                                 "-")
    
    #Ensembl Feature.query
    data$Feature.query <- ifelse(!is.na(data$Feature.query), 
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",data$Gene.query,";t=",data$Feature.query), "' target='_blank'>", data$Feature.query, "</a>")),
                                 "-")
    
    #Ensembl Feature.homolog
    data$Feature.homolog<- ifelse(!is.na(data$Feature.homolog), 
                               (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",data$Gene.homolog,";t=",data$Feature.homolog), "' target='_blank'>", data$Feature.homolog, "</a>")),
                               "-")
    
    data <- cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', data )
    
    
    
  } else if (all(c("ID.paralog","Feature.paralog" ,"SYMBOL.paralog") %in% names(data))) {
    
    
    #ClinVarID ID.query
    data$ID.query <- ifelse(data$ID.query!="NA",
                            (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",data$ID.query,"/"), "' target='_blank'>", data$ID.query, "</a>")),
                            "-")
    #ClinVarID ID.paralog
    data$ID.paralog <- ifelse(data$ID.paralog!="NA",
                              (paste0("<a href='", paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/",data$ID.paralog,"/"), "' target='_blank'>", data$ID.paralog, "</a>")),
                              "-")
    
    #Ensembl SYMBOL.query
    data$SYMBOL.query <- ifelse(!is.na(data$SYMBOL.query), 
                                (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",data$Gene.query), "' target='_blank'>", data$SYMBOL.query, "</a>")),
                                "-")
    
    #Ensembl SYMBOL.paralog
    data$SYMBOL.paralog <- ifelse(!is.na(data$SYMBOL.paralog), 
                                  (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",data$Gene.paralog), "' target='_blank'>", data$SYMBOL.paralog, "</a>")),
                                  "-")
    
    #Ensembl Feature.query
    data$Feature.query <- ifelse(!is.na(data$Feature.query), 
                                 (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",data$Gene.query,";t=",data$Feature.query), "' target='_blank'>", data$Feature.query, "</a>")),
                                 "-")
    
    #Ensembl Feature.paralog
    data$Feature.paralog <- ifelse(!is.na(data$Feature.paralog),
                                   (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=",data$Gene.paralog,";t=",data$Feature.paralog), "' target='_blank'>", data$Feature.paralog, "</a>")),
                                   "-")
    
    #print(paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",map[unlist(paralog$Gene.query)],";g1=",map[unlist(paralog$SYMBOL.paralog)]))
    #Ensembl alignment URL
    # https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=ENSG00000213281;g1=ENSG00000133703;seq=cDNA
    data$Ensembl_alignment_link <- ifelse(!is.na(data$Gene.paralog), 
                                          (paste0("<a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog/Alignment?db=core;g=",data$Gene.query,";g1=",data$Gene.paralog),"' class='btn btn-default btn-xs btn-block active' style='border-radius:1px;padding:0px' target='_blank'>alignment</a>")),
                                          "-") 
    
    data <- cbind(' ' = '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css"> <i class="fa fa-plus-square fa-lg"></i>', data )
    
    
    } else if ("Gene.paraloc" %in% names(data)) {
      
      data$Gene_ENSG.query   <- map_ENSG[unlist(data$Gene.query)] 
      data$Gene_ENSG.paraloc <- map_ENSG[unlist(data$Gene.paraloc)]
      
      #Ensembl Gene.query
      data$Gene.query <- ifelse(data$Gene.query!="NA",
                                (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",data$Gene_ENSG.query), "' target='_blank'>", data$Gene.query, "</a>")),
                                "-")
      #Ensembl Gene.paraloc
      data$Gene.paraloc <- ifelse(data$Gene.paraloc!="NA",
                                  (paste0("<a href='", paste0("https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",data$Gene_ENSG.paraloc), "' target='_blank'>", data$Gene.paraloc, "</a>")),
                                  "-")
      
      data <- select(data,-c("Gene_ENSG.query","Gene_ENSG.paraloc"))
    }
  
  return(data)
  
}






# rename DT_colnames ----

# paralog
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
                         'cdot.query' = 'cdot.query',
                         'HGVSc.query' = 'HGVSc.query',
                         'HGVSp.query' = 'HGVSp.query',
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
                         'cdot' = 'cdot.paralog',
                         'HGVSc' = 'HGVSc.paralog',
                         'HGVSp' = 'HGVSp.paralog',
                         'Residue' = 'Amino_acids.paralog',
                         'cDNA position' = 'cDNA_position.paralog',
                         'Residue position' = 'AA_position.paralog',
                         'Codons' = 'Codons.paralog',
                         'Para_Z Score' = 'Para_z_score.paralog',
                         'Ensembl alignment' = 'Ensembl_alignment_link')


# paraloc
paraloc_DT_colnames <-  c('POS' = 'CHR.query', 
                          'CHR' = 'POS.query', 
                          'REF' = 'REF.query', 
                          'Query variant' = 'var_id.query',
                          'Query gene' = 'Gene.query',
                          'Query residue' = 'AA.query',
                          'Paralogous gene' = 'Gene.paraloc',
                          'Chromosome' = 'CHR.paraloc',
                          'AA positions' = 'AA_pos.paraloc')

# homolog
homolog_DT_colnames <- c('Chr.query' = 'CHR.query',
                         'Pos.query' = 'POS.query',
                         'REF.query' = 'REF.query',
                         'ALT.query' = 'ALT.query',
                         'Query variant' = 'var_id.query',
                         'ClinVar.query' = 'ID.query',
                         'ClinVar Class.query' = 'CLNSIG.query',
                         'Gene.query' = 'SYMBOL.query',
                         'ENSG.query' = 'Gene.query',
                         'ENST.query' = 'Feature.query',
                         'cdot.query' = 'cdot.query',
                         'HGVSc.query' = 'HGVSc.query',
                         'HGVSp.query' = 'HGVSp.query',
                         'Residue.query' = 'Amino_acids.query',
                         'cDNA position.query' = 'cDNA_position.query',
                         'Residue_position.query' = 'AA_position.query',
                         'Codons.query' = 'Codons.query',
                         'Pfam domain' = 'Pfam_domain.query',
                         'Pfam position' = 'Pfam_pos.query',
                         'Chr' = 'CHR.homolog',
                         'Pos' = 'POS.homolog',
                         'REF' = 'REF.homolog',
                         'ALT' = 'ALT.homolog',
                         'Homologous variant' = 'var_id.homolog',
                         'ClinVar ID' = 'ID.homolog',
                         'ClinVar Class' = 'CLNSIG.homolog',
                         'Gene' = 'SYMBOL.homolog',
                         'ENSG' = 'Gene.homolog',
                         'ENST' = 'Feature.homolog',
                         'cdot' = 'cdot.homolog',
                         'HGVSc' = 'HGVSc.homolog',
                         'HGVSp' = 'HGVSp.homolog',
                         'Residue' = 'Amino_acids.homolog',
                         'cDNA position' = 'cDNA_position.homolog',
                         'Residue position' = 'AA_position.homolog',
                         'Codons' = 'Codons.homolog')


# childrow_JS_callback ----

childrow_JS_callback_paralog <- c("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<table style=\"border-spacing:50px\" style=\"width:50%\" cellpadding=\"50px\" style=\"padding-left:50px\">'+
      '<thead>'+
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
            '<th>ENST</th>'+
            '<td>'+ d[10] +'</td>'+
        '</tr>'+
        '<tr>'+        
            '<th>HGVSc</th>'+
    //      '<td>'+ d[9] + '(' + d[8] + '):' + d[10] + ' (' + d[11] + ')' +'</td>'+
    //      '<td>'+ d[10] + ':' + d[11] + ' (' + d[12] + ')' +'</td>'+
            '<td>'+ d[12] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>HGVSp</th>'+
            '<td>'+ d[13] +'</td>'+
        '</tr>'+
        '<tr>'+        
            '<th>Codons</th>'+
            '<td>'+ d[17] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Para_Z Score</th>'+
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


childrow_JS_callback_homolog <- c("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<table style=\"border-spacing:50px\" style=\"width:50%\" cellpadding=\"50px\" style=\"padding-left:50px\">'+
      '<thead>'+
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
            '<th>ENST</th>'+
            '<td>'+ d[10] +'</td>'+
        '</tr>'+
        '<tr>'+ 
            '<th>HGVSc</th>'+
            '<td>'+ d[12] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>HGVSp</th>'+
            '<td>'+ d[13] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Codons</th>'+
            '<td>'+ d[17] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Pfam domain</th>'+
            '<td>'+ d[18] +'</td>'+
        '</tr>'+
        '<tr>'+
            '<th>Pfam domain position</th>'+
            '<td>'+ d[19] +'</td>'+
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
    list(visible = FALSE, targets = c(1:4,6:22,27,29,32:34)),
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
    list(visible = FALSE, targets = c(0:2)),
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
    list(visible = FALSE, targets = c(1:4,6:23,28,30,33:35)),
    list(width = "100px",targets = 5),
    list(orderable = FALSE, className = 'details-control', targets = 0)))

#######
