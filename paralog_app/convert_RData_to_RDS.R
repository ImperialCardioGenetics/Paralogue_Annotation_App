for (i in c(1:22,"X","Y")){
  load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
  saveRDS(Total_annotations, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Total_annotations_chrom_",i,"_noQC.RDS"), version = 2)
  load(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RData")) #load in paralogous variant data
  saveRDS(Paraloc, file = paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/data/chrom_",i,"/Para_locations_chrom_",i,"_noQC.RDS"), version = 2)
}