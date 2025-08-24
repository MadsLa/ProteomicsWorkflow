source("source.R")
source("functions.R")
file = "raw/20220909_144243_20210817_XP1_YABH_MouseHF_PlasmaProt_Top2MDe_ProtPepMDBL_Report.tsv"

dfBGS <- data.table::fread(file, nrows=1)

MDBL_CheckFile <- function(file = "raw/20220909_144243_20210817_XP1_YABH_MouseHF_PlasmaProt_Top2MDe_ProtPepMDBL_Report.tsv"){
  positiveList <- c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.AllOccurringOrganisms")
  
  data.table::fread(file, nrows=1) %>% 
    names() %>% 
    as_tibble() -> dfTemp
  as_tibble(c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.AllOccurringOrganisms", "Negative_TEST")) %>% 
    rename("Needed in dataset" = value) %>% 
    rowwise() %>% 
    mutate("Found in dataset" = as.logical(sum(str_detect(dfTemp$value, `Needed in dataset`)))) %>% 
    gt()
}


MDBL_CheckFile()
