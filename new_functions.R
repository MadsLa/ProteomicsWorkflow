source("source.R")

importFile <- "raw/20250826_082431_20250620_YABH_JAK1_RPTEC_Exp2-Cocktail-ASO_cells_Report.tsv"
importCandidates <- "raw/20250826_081536_20220329_QE5_YABH_PlasmaProt_Groningen_pEF-Healthy_SP_MQ+hybrid search_Report_candidates.tsv"

MDBL_CheckBGSReport <- function(file){
  positiveList <- c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.StrippedSequence", "EG.ModifiedSequence", "PEP.AllOccurringOrganisms", "FALSE_TEST")
  
  dfTemp <- data.table::fread(input = file, nrows=1, showProgress=T) %>% 
    names() %>% 
    as_tibble() 
  
  as_tibble(positiveList) %>% 
    rename("Needed in dataset" = value) %>% 
    rowwise() %>% 
    mutate("Found in dataset" = as.logical(sum(str_detect(dfTemp$value, `Needed in dataset`))))-> p
  p %>% gt() -> q 
  print(q)

  if((nrow(p) - sum(p$`Found in dataset`))==1){
    print("All columns found in BGS report")
  }
}
  
MDBL_CheckBGSReport(importFile)


MDBL_LoadBGSReport <- function(file, n=Inf){
  dfTemp <- data.table::fread(input = file, nrows = n, showProgress = T) %>% 
    clean_names()
  
  data.frame("File" = file,
             "FileDate" = file.info(file)$mtime,
             "SizeGB" = paste(round(file.info(file)$size/1024^3,3), " GB"),
             "FASTA" = dfTemp %>% distinct(pg_fasta_name) %>% separate_rows(pg_fasta_name, sep =";") %>% distinct(pg_fasta_name) %>% summarise(pg_fasta_name = paste(pg_fasta_name, collapse = ", ")),
             "Species" = dfTemp %>% distinct(pep_all_occurring_organisms) %>% summarise(pep_all_occurring_organisms = paste(pep_all_occurring_organisms, collapse = "; ")),
             "Errors" = dfTemp %>% head(1) %>% pull(e_errors),
             "Warnings" = dfTemp %>% head(1) %>% pull(e_warnings),
             "Raw files" = dfTemp %>% distinct(r_file_name) %>% nrow(),
             "Conditions" = dfTemp %>% distinct(r_condition) %>% summarise(r_condition = paste(r_condition, collapse = "; "))) -> r
  
  print(as.data.frame(t(r)) %>% 
          rownames_to_column("Parameter") %>% 
          rename("Value" = "V1") %>%
          mutate(Parameter = ff(Parameter, c("pg_fasta_name", "pep_all_occurring_organisms", "Raw.files", "r_condition"), c("FASTA", "Species", "Raw files", "Conditions"), Parameter)) %>% 
          gt())  
  
  return(dfTemp)
}

MDBL_LoadBGSReport(importFile, n=100) 


MDBL_CheckCandidatesReport <- function(file){
  positiveList <- c("Genes", "ProteinGroups", "Pvalue", "Qvalue", "AVG Log2 Ratio", "FALSE_Check")
  
  dfTemp <- data.table::fread(input = file, nrows=1, showProgress=T) %>% 
    names() %>% 
    as_tibble() 
  
  as_tibble(positiveList) %>% 
    rename("Needed in dataset" = value) %>% 
    rowwise() %>% 
    mutate("Found in dataset" = as.logical(sum(str_detect(dfTemp$value, `Needed in dataset`))))-> p
  p %>% gt() -> q 
  print(q)
  
  if((nrow(p) - sum(p$`Found in dataset`))==1){
    print("All columns found in candidates file")
  }
}

MDBL_CheckCandidatesReport(importCandidates)


MDBL_LoadCandidateesReport <- function(file, n=Inf){
  dfTemp <- data.table::fread(input = file, nrows = n, showProgress = T) %>% 
    clean_names()
  
  r <- data.frame("File" = file,
             "FileDate" = file.info(file)$mtime,
             "SizeGB" = paste(round(file.info(file)$size/1024^3,3), " GB"),
             "Comparisons" = dfTemp %>% distinct(comparison_group1_group2) %>% separate_rows(comparison_group1_group2, sep =";") %>% distinct(comparison_group1_group2) %>% summarise(comparison_group1_group2 = paste(comparison_group1_group2, collapse = "; ")))

  print(as.data.frame(t(r)) %>% 
          rownames_to_column("Parameter") %>% 
          rename("Value" = "V1") %>%
          #mutate(Parameter = ff(Parameter, c("pg_fasta_name", "pep_all_occurring_organisms", "Raw.files", "r_condition"), c("FASTA", "Species", "Raw files", "Conditions"), Parameter)) %>% 
          gt())  
  
  dfTemp %>% 
    rename("padj" = "qvalue",
           "pvalue" = "pvalue",
           "gene_symbol" = "genes",
           "Uniprot_Accession" = "protein_groups") %>% 
    mutate(FoldChange = 2^avg_log2_ratio) %>% 
    select(comparison_group1_group2, Uniprot_Accession, gene_symbol, padj, pvalue, FoldChange) -> p
  p <- split(p, p$comparison_group1_group2)
    
  return(p)
}

dfCandidates <- MDBL_LoadCandidateesReport(importCandidates)

