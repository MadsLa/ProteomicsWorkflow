source("source.R")
source("functions.R")
file = "ReportBGS_Test.tsv"
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



MDBL_NameShortening2 <- function(data){
  data %>% 
    distinct(R.FileName) -> t
  
  # Split strings into a list of components
  split_strings <- str_split(t$R.FileName, "_", simplify = TRUE)
  
  unique_check <- split_strings %>%
    as_tibble() %>% 
    summarise(across(everything(), ~ n_distinct(.) == 1))
  
  # Extract column names that have only one unique value
  columns_with_one_value <- names(unique_check)[unique_check == TRUE]
  
  split_strings %>% 
    as_tibble() %>% 
    select(-paste(columns_with_one_value)) -> p
  
  p %>% 
    unite("R.FileName", 1:ncol(p), remove = FALSE)
}