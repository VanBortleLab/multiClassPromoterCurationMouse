rm(list = ls())
require(data.table)
require(dplyr)

# process file from ENSEMBL
a = fread("ensembl_export.tsv",sep="\t",header = T)
colnames(a) = c("Chr", "Start", "Stop", "URSID", "Strand", "Protein_ID", "Gene_Name", "Gene_Description", "RNA_type")
a = a[, Databases := ""]
a = a[, TOTAL_INFO := do.call(paste, c(.SD, sep = "|")), .SDcols = c("Protein_ID", "Gene_Name", "Gene_Description")]
a = a[, c("Protein_ID", "Gene_Name", "Gene_Description") := NULL]

a = a %>% filter(!grepl("_",Chr))
table(a$Chr)

x = a
x = x %>% 
  mutate(start_median=Start,stop_median=Stop) %>% 
  select(Chr,Start,Stop,URSID,Strand,RNA_type,Databases,start_median,stop_median,everything())
columns = colnames(x)
db_used = colnames(x)[10:ncol(x)]
# separate strands and elongate/shorten interval to 120
plus = subset(x, Strand == "1")
plus$Stop = as.integer(plus$Start) + 120
minus = subset(x,Strand == "-1")
minus$Start = as.integer(minus$Stop)-120
for (t1 in list("plus","minus"))
  {
    t2 = eval((parse(text=paste(t1))))
    if (nrow(t2) < 1) next
    # sort, merge, intersect 
    sorted = t2[order(t2$Chr, t2$Start), , drop=FALSE] #sorting in R because bedtools showed up error while sorting
    write.table(sorted, file="sorted.bed", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
    system(" bedtools merge -i sorted.bed > merged.bed")
    system(" bedtools intersect -a merged.bed -b sorted.bed -wb -wa > intersected.bed")
    
    # concatenating information in intersected bed file 
    y = read.table(file="intersected.bed",sep='\t',quote="")
    
    # determine which code to use based on strand direction 
    if(t1=="plus") { 
      code='z=y %>% group_by(V1,V2,V3) %>% summarize(,TSS=paste(V5,collapse="$"),start_median=median(V11),stop_median=median(V12),URSID=paste(na.omit(V7),collapse="$"),Strand=V8[1],RNA_type=V9[1]'
    } else {
      code='z=y %>% group_by(V1,V2,V3) %>% summarize(,TSS=paste(V6,collapse="$"),start_median=median(V11),stop_median=median(V12),URSID=paste(na.omit(V7),collapse="$"),Strand=V8[1],RNA_type=V9[1]'
    }
    
    for(i2 in (1:length(db_used)))
    {
      code = paste0(code, ",", db_used[i2], '=paste(unique(V',12+i2,'),collapse="$")')
      if(i2 == length(db_used)) code=paste0(code,")")
    }
    
    eval(parse(text=code)) # executing code 

    # shorten 3' end 
    if (t1 == "plus") { z$V3=z$V3-120 }
    else { z$V2=z$V2+120 }
    
    
    ##collapsing start_median and stop_median, to TSS_median... coz start and stop are related by +-120. TSS is either start or stop
    z = z %>%
      mutate(TSS_median=as.integer(ifelse(t1=="plus",start_median,stop_median))) %>%
      mutate(STOP_median=as.integer(ifelse(t1=="minus",start_median,stop_median))) %>%
      select(V1,V2,V3,TSS,TSS_median,STOP_median,everything()) %>%
      select(-c(start_median,stop_median))
    
    
    write.table(z, file = ifelse(t1 == "plus", "pfinal.bed", "mfinal.bed"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
  }
# concatenating and sorting 
system( "cat pfinal.bed mfinal.bed > cat.bed")
ff2 = read.table("cat.bed", sep="\t", header=FALSE,quote="")
colnames(ff2) = c("Chr", "Start", "Stop", "TSS", "TSS_median", "STOP_median", "URSID", "Strand", "RNA_type", db_used)
sorted = ff2 %>% 
  select(-c(STOP_median)) %>% 
  arrange(Chr, Start)  #here we remove STOP_median
if(sum(grepl("\\$",sorted$TOTAL_INFO))!=0) print("Somewhere 2 or more different PCGs are collapsed together. Try finding dollar sign in TOTAL_INFO")

# write to file
write.table(sorted, "ENSEMBL_processed.bed", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, append=F)
system("rm -r mfinal.bed pfinal.bed cat.bed merged.bed sorted.bed intersected.bed") # removing intermediate files



