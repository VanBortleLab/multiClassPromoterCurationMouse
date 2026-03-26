library(here)
library(data.table)
library(dplyr)

setwd(here::here())
source("./Step1_RNACENTRAL_column_information.R")

# PART A: create filter lists 
# load in and reformat mouse bed file 
g38_raw <- read.table("mus_musculus.GRCm39.bed", sep = "\t", quote = "", comment.char = "")
g38 <- g38_raw[, c(1, 2, 3, 4, 6, 14, 15)]
colnames(g38) <- c("Chr", "Start", "Stop", "URSID", "Strand", "RNA_type", "Databases")
g38 <- g38 %>% mutate(RNA_type = tolower(RNA_type))
g38 <- g38 %>% mutate(Databases = tolower(Databases))
g38 <- arrange(g38, Chr, Start)

# store unique RNA types 
rna_types <- unique(g38$RNA_type)
all_db    <- unique(unlist(strsplit(g38$Databases, ",")))

# split g38 and database lists by RNA type 
g38_by_type <- split(g38, g38$RNA_type)

db_by_type <- lapply(rna_types, function(type) {
  unique(unlist(strsplit(g38$Databases[g38$RNA_type == type], ",")))
})

names(db_by_type) <- rna_types
print(all_db)

# check that db files and imp_ variables from source exist

missing_files <- all_db[!file.exists(file.path("database", paste0(all_db, ".tsv")))]
if (length(missing_files) > 0)
  warning("Missing database files: ", paste(missing_files, collapse = ", "))

missing_vars <- all_db[!sapply(all_db, function(db) exists(paste0("imp_", db)))]
if (length(missing_vars) > 0)
  warning("Missing imp_ variables for: ", paste(missing_vars, collapse = ", "))

# identify important databases (filter out databases without an imp variable)
imp_db   <- all_db[sapply(all_db, function(db) {
  v <- get0(paste0("imp_", db)); !is.null(v) && length(v) > 0
})]
unimp_db <- setdiff(all_db, imp_db)

message("Important databases: ",   paste(imp_db,   collapse = ", "))
message("Unimportant databases: ", paste(unimp_db, collapse = ", "))

# filter db_by_type to only include important databases
imp_db_by_type <- lapply(db_by_type, function(dbs) intersect(dbs, imp_db))

# PART B: load RNA databases 

dir.create("databases_processed", showWarnings = FALSE)

# pull columns noted in [imp_db] variable (3 cols for each database - URS ID, RNA type, and RNA info)
load_and_trim_db <- function(db_name) {
  message("Loading: ", db_name)
  path <- file.path("database", paste0(tolower(db_name), ".tsv"))
  raw  <- fread(path, sep = "\t", quote = "")
  
  imp_cols <- get(paste0("imp_", db_name))  # column indices or names from source()d file
  trimmed  <- raw[, ..imp_cols]             # data.table column selection
  setnames(trimmed, c("URSID", "RNA_type", "Info"))
  
  # Drop empty Info rows, then collapse duplicate URSIDs
  trimmed <- trimmed[Info != ""]
  trimmed <- trimmed[, .(RNA_type = RNA_type[1],
                          Info     = paste(unique(Info), collapse = ";")),
                     by = URSID]
  message("  Rows after collapsing: ", nrow(trimmed))
  trimmed
}

db_data <- setNames(lapply(imp_db, load_and_trim_db), imp_db)
db_data$RNA_type = tolower(db_data$RNA_type)

# save cleaned databases to file 
for (db_name in imp_db) {
  write.table(db_data[[db_name]],
              file      = file.path("databases_processed", db_name),
              row.names = FALSE, quote = FALSE, sep = "\t")
}

# PART C: build masterfiles 

dir.create("rnacentral_masterfiles_raw", showWarnings = FALSE)

for (type in rna_types) {
  message("Making master file for: ", type)
  master <- g38_by_type[[type]]
  
  for (db_name in imp_db_by_type[[type]]) {
    # create subset for current RNA type 
    db_subset <- db_data[[db_name]][db_data[[db_name]]$RNA_type == type, ] 
    
    if (nrow(db_subset) == 0) next # skip if there's no lines in the processed database that match 
    
    join_data <- db_subset[, c("URSID", "Info")]
    join_data$URSID <- paste0(join_data$URSID, "_10090") 
    colnames(join_data)[2] <- db_name
    
    master <- left_join(master, join_data, by = "URSID") 
  }
  
  master <- master[, colSums(!is.na(master)) > 0]
  
  write.table(master,
              file      = file.path("rnacentral_masterfiles_raw", paste0(type, "_masterfile.bed")),
              quote     = FALSE, row.names = FALSE, sep = "\t")
}
