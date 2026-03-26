rm(list = ls())
require(dplyr)
require(data.table)


  # RNACENTRAL
  rnacentral = fread("./rnacentral_masterfiles_final/TOTAL_rnacentral.bed", sep="\t")
  rnacentral = rnacentral[!grepl("_",rnacentral$Chr),] ; rnacentral=rnacentral[!grepl("\\.",rnacentral$Chr),] ## removing bad chromosome with "_" and "."
  
  # PCG
  pcg = fread("./ENSEMBL_processed.bed", sep="\t")
  
  # combining 
  
  final = rbind(rnacentral,pcg)  
  final = as.data.frame(final) %>% rowwise() %>%
    mutate(V1=Chr, 
           V2=ifelse(Strand=="+",TSS_median-50,TSS_median-150),
           V3=ifelse(Strand=="+",TSS_median+150,TSS_median+50)) %>% 
    select(V1,V2,V3,everything())
  final = final %>% 
    rowwise() %>% 
    mutate(V2=max(0,V2))  ##Making negative V2 as 0
  final$V2 = as.integer(final$V2);final$V3=as.integer(final$V3) ##Making V2 and V3 as integer
  final = final %>% 
    dplyr::rename(TSS_Start=Start,TSS_Stop=Stop) %>% 
    dplyr::arrange(V1,V2)
  
  fwrite(final,"MultiClassPromoter.tsv", row.names = F, col.names = T, sep= "\t", quote = F,append = F)
  remove(rnacentral, pcg)

  
  

