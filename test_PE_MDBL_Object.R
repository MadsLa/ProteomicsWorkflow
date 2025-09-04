source("source.R")
source("new_functions.R")

gc() # clear memory

importBGS <- "raw/20250826_081536_20220329_QE5_YABH_PlasmaProt_Groningen_pEF-Healthy_SP_MQ+hybrid search_Report.tsv"
importCandidates <- "raw/20250826_081536_20220329_QE5_YABH_PlasmaProt_Groningen_pEF-Healthy_SP_MQ+hybrid search_Report_candidates.tsv"

#---- Check files
MDBL_CheckBGSReport(importBGS)
MDBL_CheckCandidatesReport(importCandidates)

#---- Load file
dfBGS <- MDBL_LoadBGSReport(importBGS, n = Inf)
dfCandidates <- MDBL_LoadCandidateesReport(importCandidates)

#---- Proteins
dfBGS %>% 
  filter(pg_quantity>0) %>% 
  filter(!is.na(pg_genes)) %>% 
  distinct(r_file_name, pg_protein_accessions, pg_genes, pg_quantity) %>% 
  pivot_wider(names_from = r_file_name, values_from = pg_quantity) %>% 
  rename("Uniprot_Accession" = "pg_protein_accessions",
         "gene_symbol" = "pg_genes") -> proteins

#---- Peptides
dfBGS %>% 
  filter(fg_quantity>0) %>% 
  filter(!is.na(pg_genes)) %>% 
  distinct(r_file_name, pg_protein_accessions, pg_genes, pep_stripped_sequence, fg_id, fg_quantity) %>% 
  mutate(fg_id = str_remove_all(fg_id, "_"),
         fg_id = str_replace_all(fg_id, "\\[Carbamidomethyl \\(C\\)\\]", "\\(UniMod:4\\)"),
         fg_id = str_replace_all(fg_id, "\\[Acetyl \\(Protein N\\-term\\)\\]", "\\(UniMod:1\\)"),
         fg_id = str_replace_all(fg_id, "\\[Oxidation \\(M\\)\\]", "\\(UniMod:35\\)"),
         "PEP_Count" = 1,
         "Panel_7k" = 1) %>% 
  rename("UniprotAccession" = "pg_protein_accessions",
         "Sequence_NoMods" = "pep_stripped_sequence",
         "Sequence_WithMods" = "fg_id",
         "gene_symbol" = "pg_genes") %>% 
  pivot_wider(names_from = r_file_name, values_from = fg_quantity) -> peptides

#---- Meta data
dfBGS %>% 
  distinct(r_condition, r_file_name) %>% 
  rename("Condition" = "r_condition",
         "Sample_id" = "r_file_name") -> metadata

#---- DGE
dge <- MDBL_LoadCandidateesReport(importCandidates)


#---- Experiment Meta data
info_project = "Prepare_for_ProteoExplorer_DummyData"
Project_tag = "TestFilename"
Labels = ""
MS_method = "LCMS proteomics"
Experimental_description = "Prepare_for_ProteoExplorer_DummyData"
Project = "Other_137"
Tissue = "Blood Plasma"
BRENDA = "BTO:0000131" #BRENDA Tissue Ontology https://www.ebi.ac.uk/ols/ontologies/bto
Organism = "Mouse"
Taxonomy = "10090" #https://www.uniprot.org/taxonomy/10090

info_dataset<-as_tibble(data.frame(Labels,MS_method,Experimental_description,Project,Tissue,BRENDA,Organism,Taxonomy))

info_dataset$Labels = paste0(info_dataset$Labels,Project_tag,".Rdata")


#---- Build object

# All data
data_prot1 = list()
data_prot1[["exp_data"]] = list()
data_prot1$exp_data[["normalized_intensities"]] = proteins
data_prot1[["peptides"]] = peptides
data_prot1[["meta_data"]] = metadata
data_prot1[["dge_data"]] = dge
data_prot1[["info_dataset"]] = info_project

save(data_prot1, file = paste0("PE_Output/",info_dataset$Labels))
write_delim(info_dataset,file = paste0("PE_Output/20220624MDBL_info_datasets.txt"), append = TRUE, quote = "all", delim = "\t")
