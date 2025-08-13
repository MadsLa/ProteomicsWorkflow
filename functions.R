source("source.R")

# SummaryTableSamples
'Makes simple count per filename'
SummaryTableSamples <- function(data, exportTable = F){
  data %>% 
    distinct(R.FileName, R.Condition, PG.ProteinAccessions, PG.Quantity) %>% 
    group_by(R.FileName) %>% 
    summarise("SummedIntensity" = sum(PG.Quantity, na.rm = T)) -> temp
  
  data %>% 
    group_by(R.FileName, R.Condition) %>% 
    summarise("UniqueProteins" = length(unique(PG.ProteinAccessions)),
              "UniquePeptides" = length(unique(PEP.StrippedSequence)),
              "MeanMissedCleavage" = mean(PEP.NrOfMissedCleavages)) %>% 
    left_join(temp, by = join_by(R.FileName)) %>% 
    ungroup()%>% 
    arrange(desc(R.Condition)) %>% 
    mutate(n=row_number(),
           R.FileName = fct_reorder(R.FileName, n)) %>% 
    select(-n) -> n
  
  if(exportTable==T){
    write_csv(n, file="SummarySamples.csv")
  }
  
  return(n)
}


# SummaryTableConditions
'Make simple table on counts based on conditions'
SummaryTableConditions <- function(data){
  data %>% 
    group_by(R.FileName, R.Condition) %>% 
    summarise("UniqueProteins" = length(unique(PG.ProteinAccessions)),
              "UniquePeptides" = length(unique(PEP.StrippedSequence))) %>% 
    group_by(R.Condition) %>% 
    summarise("Files" = n(), 
              "MeanProteins" = mean(UniqueProteins), 
              "RangeProteins" = paste0(min(UniqueProteins) , "-", max(UniqueProteins)),
              "MeanPeptides" = mean(UniquePeptides), 
              "RangePeptides" = paste0(min(UniquePeptides) , "-", max(UniquePeptides))) %>% 
    rename("Condition" = R.Condition)
}


MDBL_HeatMap <- function(data){
  data %>% 
    select(R.FileName, R.Condition, PG.ProteinAccessions, PG.Quantity) %>% 
    distinct() %>% 
    mutate(PG.Quantity = log(PG.Quantity)) %>% 
    filter(PG.Quantity>1) %>% 
    pivot_wider(id_cols  = c(PG.ProteinAccessions), names_from = R.FileName, values_from = PG.Quantity) %>% 
    column_to_rownames("PG.ProteinAccessions") %>% 
    as.data.frame() -> dfTemp
  
  dfTemp %>% 
    cor(method = "pearson", use="pairwise.complete.obs") -> corr
  
  data.frame(data %>% select(R.FileName, R.Condition) %>% distinct()) %>% column_to_rownames("R.FileName") -> annot
  
  pheatmap::pheatmap(corr,
                     show_rownames = F, 
                     annotation_col = annot,
                     show_colnames = F, 
                     fontsize_row = 1,
                     fontsize_col = 6,
                     width = 8)
  
}

MDBL_PlotCV <- function(data){
  data %>% 
    distinct(R.FileName, R.Condition, PG.ProteinAccessions, PG.Quantity) %>% 
    #filter(PG.Quantity>1) %>% 
    mutate(PG.Quantity = log(base=10, PG.Quantity)) %>% 
    group_by(R.Condition, PG.ProteinAccessions) %>% 
    summarise(n=n(),
              cv = CV(PG.Quantity)) %>% 
    filter(cv>0) %>% 
    ggplot(aes(cv, R.Condition, fill=R.Condition)) + 
    geom_violin(draw_quantiles = T)+
    scale_x_continuous(labels = scales::percent)+
    labs(title="CV on Log10 Intensities", x="CV", y="")+
    theme_minimal()+
    theme(legend.position = "none")
}

MDBL_ProteinsPerFileBarplot <- function(data){
  data %>% 
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
}

MDBL_PlotIntensityDensity <- function(data){
  data %>% 
    select(R.FileName, R.Condition, PG.ProteinAccessions, PG.Quantity) %>% 
    distinct() %>% 
    filter(PG.Quantity>1) %>% 
    mutate(PG.Quantity = log(base=2, PG.Quantity)) %>% 
    ggplot(aes(PG.Quantity, group = R.FileName, color=R.Condition))+
    geom_density()+
    facet_wrap(~R.Condition)+
    labs(title="Log10 Quantity Distribution", x="", y="log10(Quantity")+
    theme_minimal()+
    theme(legend.position = "none")
}

MDBL_PlotMissedCleavage <- function(data){
  data %>% 
    mutate(PEP.NrOfMissedCleavages = fct_rev(as.factor(PEP.NrOfMissedCleavages))) %>% 
    group_by(R.FileName, R.Condition, PEP.NrOfMissedCleavages) %>% 
    tally() %>% 
    ggplot(aes(n, R.FileName, fill=PEP.NrOfMissedCleavages))+
    geom_col(position = "fill", width=1)+
    scale_x_continuous(labels = scales::percent)+
    scale_fill_brewer(name = "Missed \ncleavages", palette = "Spectral")+
    facet_wrap(~R.Condition, scales="free_y")+
    labs(title="Missed cleavage sites per peptide", x="", y="")+
    theme_minimal()+
    theme(axis.text.y = element_blank())
}

MDBL_SamplesPerCondition <- function(data){
  data %>% 
    distinct(R.FileName, R.Condition) %>% 
    group_by(R.Condition) %>% 
    tally() %>% 
    ggplot(aes(n, R.Condition, fill=R.Condition, label=n))+
    geom_col()+
    geom_text(hjust=1.3)+
    labs(title="Samples per condition", x="samples per condition", y="")+
    theme_minimal()+
    theme(legend.position = "none")
}

MDBL_ProteinsPerFileBoxplot <- function(data){
  data %>% 
    distinct(R.Condition, R.FileName, PG.ProteinAccessions) %>% 
    group_by(R.FileName, R.Condition) %>% 
    tally() -> n
  get_box_stats <- function(y, upper_limit = max(n %>% pull(n)) * 1.15) {
    return(data.frame(
      y = 0.95 * upper_limit,
      label = paste("CV =", round(CV(y)*100, 2), "%")
    ))
  }
  n %>%
    ggplot(aes(n, R.Condition, fill=R.Condition))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(alpha=0.3, height = 0.2, pch=16)+
    stat_summary(fun.data = get_box_stats, geom = "text", vjust = 0.5, hjust = 0.9, size=3)+
    labs(title="Unique proteins per raw file", x="Unique proteins", y="")+
    theme_minimal()+
    theme(legend.position = "none",
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())
}

MDBL_FilterVV <- function(data, pct=0.3){
  pctFilter=pct
  data %>% 
    distinct(R.FileName, R.Condition) %>% 
    group_by(R.Condition) %>% 
    tally(name = "n_all") -> dfTemp
  
  data %>% 
    distinct(R.FileName, R.Condition, PG.ProteinAccessions, PG.Quantity) %>%
    group_by(R.Condition, PG.ProteinAccessions) %>%
    mutate(n_condition = n()) %>% 
    left_join(dfTemp, by = join_by(R.Condition)) %>% 
    mutate(pct = n_condition/n_all) %>% 
    filter(pct>=pctFilter) -> data
  return(data)
}

MDBL_DataCompleteness <- function(data){
  data %>% 
    distinct(R.FileName, R.Condition) %>% 
    group_by(R.Condition) %>% 
    tally(name = "SamplesTotal") -> dfTemp
  
  data %>% 
    distinct(R.Condition, PG.ProteinAccessions, R.FileName) %>% 
    group_by(R.Condition, PG.ProteinAccessions) %>% 
    tally(name ="SamplesID") %>% 
    ungroup() %>% 
    left_join(dfTemp, by = join_by(R.Condition)) %>% 
    mutate("pct" = SamplesID/SamplesTotal) %>% 
    group_by(R.Condition) %>% 
    arrange(desc(pct)) %>% 
    mutate(n=row_number())  %>% 
    ungroup() %>%
    ggplot(aes(n, pct, color=R.Condition)) +
    geom_line(linewidth=0.75)+
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1))+
    facet_wrap(~R.Condition)+
    labs(title="Data completeness", x="proteins", y="% completeness")+
    theme_minimal()+
    theme(legend.position = "none",
          panel.grid.minor.y = element_blank())
  
}
