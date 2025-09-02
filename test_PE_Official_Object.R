## Create object to upload to proteomics app:

#### prepare for Proteomics app #################################################
### This script prepares data in the format compatible with ProteoExplorer. Modified from **KRMB** script for Somalogic proteomics data analysis
### Input data: Statistical analysis results, raw proteomics data, and experiment metadata.
### Output data: ProteoExplorer object
### Target audience: Collaborators outside of DSI who performed data analysis but lack the integrated pipeline to upload results to ProteoExplorer

library(dplyr)
library(tidyr)
library(tibble)
library(readr)

#Storage path of all ProteoExplorer objects 
output_dir = "/novo/omdb/pds01/PDS2956/data/"

#Description of data:
info_project = "Prepare_for_ProteoExplorer_DummyData"


#Experimental information:
Project_tag="Prepare_for_ProteoExplorer_DummyData"
Labels=""
MS_method="Somalogic"
Experimental_description="Prepare_for_ProteoExplorer_DummyData"
Project="Other_137"
Tissue="Blood Plasma"
BRENDA="BTO:0000131" #BRENDA Tissue Ontology https://www.ebi.ac.uk/ols/ontologies/bto
Organism="Mouse"
Taxonomy="10090" #https://www.uniprot.org/taxonomy/10090

info_dataset<-as_tibble(data.frame(Labels,MS_method,Experimental_description,Project,Tissue,BRENDA,Organism,Taxonomy))

info_dataset$Labels = paste0(info_dataset$Labels,Project_tag,".Rdata")

# Prepare metadata
# Below a dummy metadata (replace this with your experimental setup)
metadata <- data.frame(
  Condition = c("A", "B", "A", "B", "C"),
  Sample_id = c("S1", "S2", "S3", "S4", "S5")
)

# Prepare protein data matrix
# Create a dummy dataframe for proteins
proteins <- data.frame(
  Uniprot_Accession = c("P111", "P222", "P333", "P444", "P444"),
  gene_symbol = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneD"),
  SeqId = c("seq.111", "seq.222", "seq.333", "seq.444", "seq.555"),
  S1 = c(1.2, 3.4, 2.5, 1.8, 2.0),  # Replace with actual intensity values
  S2 = c(2.3, 1.5, 3.2, 2.7, 1.9),  # Replace with actual intensity values
  S3 = c(1.8, 2.1, 1.9, 3.0, 2.4),  # Replace with actual intensity values
  S4 = c(2.7, 1.9, 2.5, 1.6, 3.1),  # Replace with actual intensity values
  S5 = c(2.0, 1.6, 2.9, 1.7, 2.3)   # Replace with actual intensity values
)

#Drop NA genes names
proteins<-proteins[!is.na(proteins$gene_symbol),]
#Create new gene names for duplicates
dups<-proteins$`gene_symbol`[duplicated(proteins$`gene_symbol`)]
for(x in dups){
  proteins[proteins$`gene_symbol`==x,"gene_symbol"]<-paste(proteins[proteins$`gene_symbol`==x,"gene_symbol"],proteins[proteins$`gene_symbol`==x,"SeqId"],sep="_")
}       

#proteins_totest<-proteins[,-1]
if (!all(names(proteins)[4:length(names(proteins))] == metadata$Sample_id)){
  stop("sample names are in wrong order between proteins and annotations")
}

proteins<-as_tibble(proteins)

# DGE data
#Column names Uniprot_Accession gene_symbol FoldChange  padj

#Create a dummy statistical result
#Note ProteoExplorer converts FoldChange into log2 scale, so the provided value should not be log2-transformed 
stat_model1 <- data.frame(
  Uniprot_Accession = c("P111", "P222", "P333", "P444", "P444"),
  gene_symbol = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneD"),
  SeqId = c("seq.111", "seq.222", "seq.333", "seq.444", "seq.555"),
  FoldChange = c(1.5, 4.0, 0.5, 1.2, 1.9),  # Replace with actual fold change values
  padj = c(0.01, 0.005, 0.03, 0.02, 0.008),  # Replace with actual adjusted p-values
  p_value = c(0.002, 0.001, 0.006, 0.004, 0.0016)
)

dge_data = list()
result_df_lists_all<-list("GroupA_vs_GroupB"=stat_model1
)

align_results_with_app<-function(x){
  #Drop NA genes names
  y<-x[!is.na(x$gene_symbol),]
  #Create new gene names for dublicates
  dups<-y$`gene_symbol`[duplicated(y$`gene_symbol`)]
  for(dup in dups){
    y[y$`gene_symbol`==dup,"gene_symbol"]<-paste(y[y$`gene_symbol`==dup,"gene_symbol"],y[y$`gene_symbol`==dup,"SeqId"],sep="_")
  }       
  y<-as_tibble(y)
  return(y)
}

for (list_name in names(result_df_lists_all)){
  result_df_lists_all[[list_name]]<-align_results_with_app(result_df_lists_all[[list_name]])
  
}

#Create synthetic peptide data since it is needed by app.
peptides<-data.frame(Sequence_NoMods=rep("ABC",dim(proteins)[1]),Sequence_WithMods=rep("ABC",dim(proteins)[1]),UniprotAccession=proteins$Uniprot_Accession,gene_symbol=proteins$gene_symbol,PEP_Count=rep(1,dim(proteins)[1]),Panel_7k=rep(1,dim(proteins)[1]))

# All data
data_prot1 = list()
data_prot1[["exp_data"]] = list()
data_prot1$exp_data[["normalized_intensities"]] = proteins
data_prot1[["peptides"]] = peptides
data_prot1[["meta_data"]] = metadata
data_prot1[["dge_data"]] = result_df_lists_all
data_prot1[["info_dataset"]] = info_project

save(data_prot1,file = paste0(output_dir,"/",info_dataset$Labels))
#Update info data
#write_delim(info_dataset,file = #paste0(output_dir,"/20220624MDBL_info_datasets.txt"),append = TRUE,quote = "all",delim = "\t")