for (i in c(1:22,"X","Y")){
  load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
  saveRDS(Total_annotations, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RDS"), version = 2)
  load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
  saveRDS(Paraloc, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RDS"), version = 2)
}

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
saveRDS(raw_data, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/Total_annotations_all_chrom_noQC.RDS"), version = 2)

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
saveRDS(Paraloc_data, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/Para_locations_all_chrom_noQC.RDS"), version = 2)

test = readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/Total_annotations_all_chrom_noQC.RDS"))
test = readRDS(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/Para_locations_all_chrom_noQC.RDS"))
