source("source.R")
source("functions.R")
file = "raw/20250702_100140_20250610_YABH_MFGE8_Mouse_PlasmaProt_MagNet_PE_Report.tsv"
dfBGS <- data.table::fread(file)



dfBGS %>% 
  MDBL_tsne()
