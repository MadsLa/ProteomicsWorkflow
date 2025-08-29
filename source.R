pacman::p_load(tidyverse,
               extrafont,
               ggtext,
               showtext,
               patchwork,
               gt,
               janitor,
               ggforce)

Sys.setlocale("LC_TIME", "English")
options(scipen=999)
options(dplyr.summarise.inform = FALSE)

#list <- c("#A9AF7E", "#557153", "#D9DCCA", "#394D3F", "#CABBA2", "#7E9B76", "#E0E5C1", "#6C7B6B", "#9C9C5E", "#8A9483")


ff = function(x, patterns, replacements = patterns, fill = NA, ...)
  #find/extract function
  #use: ff(df$Rawfile, c("wt", "dbdb"), c("WT", "DBDB"), "Other", ignore.case = TRUE)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  return(ans)
}

show_in_excel <- function(.data){
  #usage: mtcars %>%  show_in_excel() 
  if(interactive()) {# avoid unwanted Excel's executions
    tmp <- tempfile(fileext = ".xlsx")
    writexl::write_xlsx(.data, tmp)
    #readr::write_excel(.data, tmp)
    fs::file_show(tmp)
  }
  
  .data
  
}

head_tail <- function(d, m=5){
  library(dplyr)
  # print the head and tail together
  HEAD = head(d,m)
  TAIL = tail(d,m)
  HEAD %>% 
    rows_append(TAIL)
}

CV <- function(x){
  (sd(x)/mean(x))
}

