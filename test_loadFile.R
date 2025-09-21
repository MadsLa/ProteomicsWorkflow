source("source.R")
source("functions.R")
file = "raw/20250826_082431_20250620_YABH_JAK1_RPTEC_Exp2-Cocktail-ASO_cells_Report.tsv"

dfBGS <- data.table::fread(file, nrows=1)

MDBL_CheckFile <- function(file = NA){
  positiveList <- c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.StrippedSequence", "EG.ModifiedSequence", "PEP.AllOccurringOrganisms", "FALSE_TEST")
  
  data.table::fread(file, nrows=1) %>% 
    names() %>% 
    as_tibble() -> dfTemp
  as_tibble(positiveList) %>% 
    rename("Needed in dataset" = value) %>% 
    rowwise() %>% 
    mutate("Found in dataset" = as.logical(sum(str_detect(dfTemp$value, `Needed in dataset`))))
}


MDBL_CheckFile(file = "raw/20250826_082431_20250620_YABH_JAK1_RPTEC_Exp2-Cocktail-ASO_cells_Report.tsv")



MDBL_ReadBGSReport <- function(file, RemoveNA = T){
  positiveList <- c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.StrippedSequence", "EG.ModifiedSequence", "PEP.AllOccurringOrganisms", "FALSE_TEST")
  
  data.table::fread(file, nrows=1, showProgress=T) %>% 
    names() %>% 
    as_tibble() -> dfTemp
  
  as_tibble(positiveList) %>% 
    rename("Needed in dataset" = value) %>% 
    rowwise() %>% 
    mutate("Found in dataset" = as.logical(sum(str_detect(dfTemp$value, `Needed in dataset`))))-> p
    p %>% gt() -> q 
  print(q)
  rm(dfTemp)
  
  if((nrow(p) - sum(p$`Found in dataset`))==1){
    print("All columns found, data being loaded")
    dfTemp <- data.table::fread(file = file, showProgress=F)
  print(paste0(length(unique(dfTemp$PG.ProteinAccessions)), " unique proteins identified across ", length(unique(dfTemp$R.FileName)), " samples in ", length(unique(dfTemp$R.Condition)) , " conditions"))
  nProteins <- length(dfTemp %>%  distinct(PG.ProteinAccessions, R.FileName, PG.Quantity) %>% filter(is.na(PG.Quantity)))
  nPeptides <- length(dfTemp %>%  filter(is.na(PG.Quantity)))
  print(paste0(nPeptides, " peptides from ", nProteins, " proteins found with protein intensity as 'NaN'"))
  if(RemoveNA == T){
    dfTemp %>% filter(!is.na(PG.Quantity)) -> dfTemp
    print("Peptides/Proteins with 'NaN' removed")
  } else {
    print("Peptides/Proteins with 'NaN' NOT removed")
  }
  
  data.frame("File" = file,
             "Size, GB" = file.info(file)$size/1024^3,
             "Date" = file.info(file)$mtime,
             "FASTA" = dfTemp %>% distinct(E.LFQMethod, PG.FASTAName) %>% group_by(PG.FASTAName) %>%  summarise(n = paste(PG.FASTAName, collapse = ";")) %>% pull(n),
             "Species" = dfTemp %>% distinct(E.LFQMethod, PEP.AllOccurringOrganisms) %>% group_by(E.LFQMethod) %>%  summarise(n = paste(PEP.AllOccurringOrganisms, collapse = ";")) %>% pull(n),
             "Errors" = dfTemp %>% head(1) %>% pull(E.Errors),
             "Warnings" = dfTemp %>% head(1) %>% pull(E.Warnings),
             "Files" = dfTemp %>% distinct(R.FileName) %>% nrow(),
             "Conditions" = dfTemp %>% distinct(R.Condition) %>% nrow()) -> r
  
  print(r %>% gt())
  return(dfTemp)
    
  }
  
}

dfBGS <- MDBL_ReadBGSReport(file=file)