source("source.R")
source("functions.R")
file = "raw/20220909_144243_20210817_XP1_YABH_MouseHF_PlasmaProt_Top2MDe_ProtPepMDBL_Report.tsv"
dfBGS <- data.table::fread(file)

lookup <- "Ttc37"


dfBGS %>% 
  filter(str_detect(PG.Genes, lookup)) %>% 
  distinct(PG.Genes, PG.ProteinAccessions, EG.StrippedSequence, EG.ModifiedSequence) -> p

if(nrow(p)>0){
  p
} else {
  dfBGS %>% 
    filter(str_detect(PG.ProteinAccessions, lookup)) %>% 
    distinct(PG.Genes, PG.ProteinAccessions, EG.StrippedSequence, EG.ModifiedSequence)
}



MDBL_Lookup <- function(data, lookup="F7AI87"){
  dfBGS %>% 
    filter(PG.Quantity>0) %>% 
    filter(PEP.Quantity>0) -> s
  s %>% 
    filter(str_detect(PG.Genes, lookup)) %>% 
    distinct(PG.Genes, PG.ProteinAccessions, EG.StrippedSequence, EG.ModifiedSequence) -> p
  
  s %>% 
    filter(str_detect(PG.ProteinAccessions, lookup)) %>% 
    distinct(PG.Genes, PG.ProteinAccessions, EG.StrippedSequence, EG.ModifiedSequence) -> q
  
  if(nrow(p)>0){
    p
  } else if(nrow(q)>0) {
    q
  } else {
    print(paste0("Lookup '", lookup, "' not found!"))
  }
  
}

MDBL_Lookup(lookup = "Igkv9-124")

lookup="Igkv9-124"
dfBGS %>% 
  filter(PEP.Quantity>0) %>% 
  filter(PG.Quantity>0) %>% 
  filter(str_detect(PG.Genes, lookup)) %>% 
  distinct(R.FileName, R.Condition, PG.Genes, EG.StrippedSequence, EG.ModifiedSequence, PEP.Quantity) %>% 
  group_by(R.FileName, R.Condition) %>% 
  distinct() %>% 
  tally() %>%
  ungroup() %>% 
  mutate(R.FileName = fct_reorder(R.FileName, R.Condition)) %>% 
  ggplot(aes(n, R.FileName, fill=R.Condition))+
  geom_col()+
  labs(title="Unique peptides per raw file", subtitle="", x="", y="")+
  theme_minimal()
