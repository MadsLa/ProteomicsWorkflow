source("source.R")
dfBGS <- data.table::fread("raw/20250620_075954_20250610_YABH_MFGE8_Mouse_PlasmaProt_MagNet_BGS_Report.tsv")
dfProtPep <- data.table::fread("raw/20250620_075954_20250610_YABH_MFGE8_Mouse_PlasmaProt_MagNet_PepProt_Report.tsv")

dfBGS %>% 
  group_by(R.FileName, R.Condition) %>% 
  tally() %>% 
  ungroup() %>% 
  ggplot(aes(R.FileName, n, fill = R.Condition))+
  geom_col(width=0.8)+
  facet_wrap(~R.Condition, scales="free_x")+
  labs(title="PSMs per raw file grouped by condition", x="", y="PSM")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_blank())

dfBGS %>% 
  distinct(R.FileName, R.Condition, PG.ProteinAccessions) %>% 
  group_by(R.FileName, R.Condition) %>% 
  tally() %>% 
  ungroup() %>% 
  ggplot(aes(R.FileName, n, fill = R.Condition))+
  geom_col(width=0.8)+
  facet_wrap(~R.Condition, scales="free_x")+
  labs(title="Protein groups per raw file grouped by condition", x="", y="PSM")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_blank())