source("source.R")
source("functions.R")
file = "raw/20250703_073723_StemCells_ANWR_Run5_Plasma_PE_Report.tsv"
dfBGS <- data.table::fread(file)


dfBGS %>% 
  distinct(R.FileName, R.Condition) -> dfTemp

dfBGS %>% 
  ungroup() %>% 
  mutate(PG.Quantity = log2(PG.Quantity)) %>% 
  distinct(R.Label, PG.ProteinAccessions, PG.Quantity) %>% 
  pivot_wider(names_from = R.Label, values_from = PG.Quantity) %>% 
  na.omit() %>% 
  t() %>% 
  row_to_names(row_number = 1) %>% 
  as_tibble(rownames = "R.FileName") %>% 
  mutate(R.FileName = str_remove(R.FileName, ".raw")) %>% 
  left_join(dfTemp, by = join_by(R.FileName)) -> p


p %>% 
  select(-c(R.FileName, R.Condition)) %>%
  mutate_if(is.character,as.numeric) -> q

pacman::p_load(Rtsne)
  
tsne_result <- Rtsne(q, dims = 2, perplexity = 20, verbose = TRUE, check_duplicates = FALSE)

tsne_data <- data.frame(tsne_result$Y)
colnames(tsne_data) <- c("Dim1", "Dim2")
tsne_data$R.Condition <- p$R.Condition  # Add species labels

ggplot(tsne_data, aes(x = Dim1, y = Dim2, color = R.Condition)) +
  geom_point(size = 3) +
  labs(title = "t-SNE Plot of Iris Dataset") +
  theme_minimal()+
  theme(legend.position = "none")
