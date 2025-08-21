source("source.R")
source("functions.R")
file = "raw/20250702_100140_20250610_YABH_MFGE8_Mouse_PlasmaProt_MagNet_PE_Report.tsv"

dfBGS <- data.table::fread(file, nrows=1)

MDBL_CheckFile <- function(file = "raw/20250702_100140_20250610_YABH_MFGE8_Mouse_PlasmaProt_MagNet_PE_Report.tsv"){
  positiveList <- c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.AllOccurringOrganisms")
  
  data.table::fread(file, nrows=1) %>% 
    names() %>% 
    as_tibble() -> dfTemp
  for (n in positiveList){
    dfTemp %>% 
      mutate(Found = str_detect(value, n))
  }
}

MDBL_CheckFile()

dfTemp %>% 
  rowwise() %>% 
  mutate("Found" = sum(str_detect(value, positiveList))) %>% 
  View()

         
as_tibble(c("PG.Genes", "PG.ProteinAccessions", "PG.ProteinNames", "E.Errors", "E.Warnings", "R.Label", "R.Condition", "PG.FASTAName", "PG.Quantity", "PEP.AllOccurringOrganisms", "Negative_TEST")) %>% 
  rename("Needed in dataset" = value) %>% 
  rowwise() %>% 
  mutate("Found in dataset" = as.logical(sum(str_detect(dfTemp$value, `Needed in dataset`)))) %>% 
  gt()

         