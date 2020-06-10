library(RSQLite)

# for (i in c(1:22,"X","Y")){
#   load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
#   saveRDS(Total_annotations, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RDS"), version = 2)
#   load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
#   saveRDS(Paraloc, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RDS"), version = 2)
# }

raw_data = NULL
for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
# for (i in c(1)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
  Total_annotations = readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RDS"))
  
  if (is.null(raw_data)){
    raw_data = Total_annotations
  } else {
    raw_data = base::rbind(raw_data, dplyr::setdiff(Total_annotations, raw_data))
  }
}
rm(Total_annotations)
Paraloc_data = NULL
for (i in c(1:22,"X","Y")){ #FOR FULL DATASET UNCOMMENT AND USE THIS LINE
# for (i in c(1)){ #FOR TEST DATASET UNCOMMENT AND USE THIS LINE
  Paraloc = readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RDS"))
  Paraloc = dplyr::distinct(Paraloc)
  if (is.null(Paraloc_data)){
    Paraloc_data = Paraloc
  } else {
    Paraloc_data = base::rbind(Paraloc_data, dplyr::setdiff(Paraloc, Paraloc_data))
  }
}
rm(Paraloc)

con <- dbConnect(RSQLite::SQLite(), paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/db.sqlite"))
dbWriteTable(con, "raw_data", raw_data)
dbWriteTable(con, "Paraloc_data", Paraloc_data)

input_data = "1 115256528 T G"
input_data2 = "1 115256528 T"
output = raw_data[raw_data$var %in%  input_data$mutation,]

subtable <- dbGetQuery(
  con, paste0("SELECT * FROM raw_data WHERE var = '",input_data,"'")
)
subtable2 = dbGetQuery(
  con, paste0("SELECT * FROM Paraloc_data WHERE var = '",input_data2,"'")
)


test = dbConnect(RSQLite::SQLite(), "")
dbWriteTable(test, "raw_data", raw_data)
dbWriteTable(test, "Paraloc_data", Paraloc_data)

paste(x$mut, collapse = '", "')

test_out = dbGetQuery(
  test, paste0("SELECT * FROM raw_data WHERE var IN ('",paste(x$mut, collapse = "','"),"')")
)
