source("source.R")
source("functions.R")
file = "raw/20240725_095500_MFGE8_BMDiscovery_directDIA_crossExperiment_LYTL_RHZZ_Report.tsv"
df <- data.table::fread(file)


MDBL_NameShortening <- function(data){
  data %>% 
  distinct(R.FileName) %>% 
  mutate(NameLength=str_length(R.FileName)) -> t

for(n in max(t$NameLength):3){
  t %>% 
    mutate(name = str_trunc(R.FileName, n, ellipsis = "")) %>%  
    distinct(name) %>% nrow() -> s
  if(s==1){
    
    t %>% 
      mutate(name = str_trunc(R.FileName, n, ellipsis = "")) %>%  
      distinct(name) %>% pull(name) -> k
    print(paste0("String to remove: ",k))
    break 
  }
  
}
  data %>% 
    mutate(R.FileName = str_remove(R.FileName, k)) -> data
  return(data)
}

df %>% 
  MDBL_NameShortening()


data <- dfBGS %>% distinct(R.FileName)

data %>% 
  mutate(parts = str_count(R.FileName, "_"),
         max_parts = max(parts, na.rm=T)) 
  
max(, na.rm = T )

    