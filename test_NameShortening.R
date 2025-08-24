source("source.R")
source("functions.R")
file = "raw/20220909_144243_20210817_XP1_YABH_MouseHF_PlasmaProt_Top2MDe_ProtPepMDBL_Report.tsv"
dfBGS <- data.table::fread(file)

dfBGS %>% 
  distinct(R.FileName) %>% 
  mutate(NameLength=str_length(R.FileName)) -> t

for(n in max(t$NameLength):3){
  dfBGS %>% 
    mutate(name = str_trunc(R.FileName, n, ellipsis = "")) %>%  
    distinct(name) %>% nrow() -> s
  if(s==1){
    print(n)
    dfBGS %>% 
      mutate(name = str_trunc(R.FileName, n, ellipsis = "")) %>%  
      distinct(name) %>% pull(name)
    break 
  }
  
}

  
