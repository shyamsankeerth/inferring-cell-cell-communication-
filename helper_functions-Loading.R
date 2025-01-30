library(shiny)
library(dplyr)
library(ggplot2)


####realdata and the gene_dictionary are the two datasets which has the entire data of all genes and interactions##############
realdata <- readRDS("data/sig_hits_all_studies.rds")
gex_list <-readRDS("~/crosstalk_webserver/data/density/gex_list.rds") #gex_list contains the geneexpression calculated for each gene per celltype and filtered by different datasets and also has cumilative for gene per each celltype
# Rename specific cell types 
realdata <- realdata %>%
  mutate(
    emitter = case_when(
      emitter %in% c("ASC1", "ASC2", "ASC3", "ASC4", "Astro", "Astro.0", "Astro.1", "Astro.2", "Astro.3", "Astro.4") ~ "  Astrocytes",
      emitter %in% c("Endo", "endo.ECs", "endo.FLCs", "endo.VSMC", "PER_END1", "PER_END2", "PER_END3") ~ "  Endothelial",
      emitter %in% c("MG1", "MG2", "MG3", "Mic", "Micro", "Micro.0", "Micro.1", "Micro.2", "Micro.3", "Micro.4", "Micro.5", "Micro.6", "Micro.7", "Micro.8") ~ "  Microglia",
      emitter %in% c("Ex1","Ex0", "EX1", "EX2", "EX3", "EX4", "EX5", "Excit", "NeuronE.0", "NeuronE.1", "NeuronE.2", "NeuronE.6") ~ "  NeuronE",
      emitter %in% c("In", "INH1", "INH2", "INH3", "INH4", "Inhit", "NeuronI.3", "NeuronI.4", "NeuronI.5") ~ "  NeuronI",
      emitter %in% c("OPC", "OPC.0", "OPC.1", "OPC.2", "OPC.3", "OPC.4", "OPC.5", "OPC.6", "OPC1", "OPC2") ~ "  OPC",
      emitter %in% c("ODC1", "ODC10", "ODC11", "ODC12", "ODC13", "ODC2", "ODC3", "ODC4", "ODC5", "ODC6", "ODC7", "ODC8", "ODC9",
                     "Oli0", "Oli1", "Oligo", "Oligo.0", "Oligo.1", "Oligo.2", "Oligo.3", "Oligo.4", "Oligo.5", "Oligo.6", "Oligo.7", "Oligo.8") ~ "  Oligos",
      emitter %in% c("U_ASC1", "U_ASC2", "U_ASC3", "U_ASC4", "U_Astro", "U_Astro.0", "U_Astro.1", "U_Astro.2", "U_Astro.3", "U_Astro.4") ~ "U_Astrocytes",
      emitter %in% c("U_Endo", "U_endo.ECs", "U_endo.FLCs", "U_endo.VSMC", "U_PER_END1", "U_PER_END2", "U_PER_END3") ~ "U_Endothelial",
      emitter %in% c("U_MG1", "U_MG2", "U_MG3", "U_Mic", "U_Micro", "U_Micro.0", "U_Micro.1", "U_Micro.2", "U_Micro.3", "U_Micro.4", "U_Micro.5", "U_Micro.6", "U_Micro.7", "U_Micro.8") ~ "U_E_Microglia",
      emitter %in% c("U_Ex1", "U_Ex0", "U_EX1", "U_EX2", "U_EX3", "U_EX4", "U_EX5", "U_Excit", "U_NeuronE.0", "U_NeuronE.1", "U_NeuronE.2", "U_NeuronE.6") ~ "U_E_NeuronE",
      emitter %in% c("U_In", "U_INH1", "U_INH2", "U_INH3", "U_INH4", "U_Inhit", "U_NeuronI.3", "U_NeuronI.4", "U_NeuronI.5") ~ "U_E_NeuronI",
      emitter %in% c("U_OPC", "U_OPC.0", "U_OPC.1", "U_OPC.2", "U_OPC.3", "U_OPC.4", "U_OPC.5", "U_OPC.6", "U_OPC1", "U_OPC2") ~ "U_E_OPC",
      emitter %in% c("U_ODC1", "U_ODC10", "U_ODC11", "U_ODC12", "U_ODC13", "U_ODC2", "U_ODC3", "U_ODC4", "U_ODC5", "U_ODC6", "U_ODC7", "U_ODC8", "U_ODC9",
                     "U_Oli0", "U_Oli1", "U_Oligo", "U_Oligo.0", "U_Oligo.1", "U_Oligo.2", "U_Oligo.3", "U_Oligo.4", "U_Oligo.5", "U_Oligo.6", "U_Oligo.7", "U_Oligo.8") ~ "U_E_Oligos",
      TRUE ~ emitter
    ),
    receiver = case_when(
      receiver %in% c("ASC1", "ASC2", "ASC3", "ASC4", "Astro", "Astro.0", "Astro.1", "Astro.2", "Astro.3", "Astro.4") ~ "Astrocytes",
      receiver %in% c("Endo", "endo.ECs", "endo.FLCs", "endo.VSMC", "PER_END1", "PER_END2", "PER_END3") ~ "Endothelial",
      receiver %in% c("MG1", "MG2", "MG3", "Mic", "Micro", "Micro.0", "Micro.1", "Micro.2", "Micro.3", "Micro.4", "Micro.5", "Micro.6", "Micro.7", "Micro.8") ~ "Microglia",
      receiver %in% c("Ex1", "Ex0","EX1", "EX2", "EX3", "EX4", "EX5", "Excit", "NeuronE.0", "NeuronE.1", "NeuronE.2", "NeuronE.6") ~ "NeuronE",
      receiver %in% c("In", "INH1", "INH2", "INH3", "INH4", "Inhit", "NeuronI.3", "NeuronI.4", "NeuronI.5") ~ "NeuronI",
      receiver %in% c("OPC", "OPC.0", "OPC.1", "OPC.2", "OPC.3", "OPC.4", "OPC.5", "OPC.6", "OPC1", "OPC2") ~ "OPC",
      receiver %in% c("ODC1", "ODC10", "ODC11", "ODC12", "ODC13", "ODC2", "ODC3", "ODC4", "ODC5", "ODC6", "ODC7", "ODC8", "ODC9",
                      "Oli0", "Oli1", "Oligo", "Oligo.0", "Oligo.1", "Oligo.2", "Oligo.3", "Oligo.4", "Oligo.5", "Oligo.6", "Oligo.7", "Oligo.8") ~ "Oligos",
      
      receiver %in% c("U_ASC1", "U_ASC2", "U_ASC3", "U_ASC4", "U_Astro", "U_Astro.0", "U_Astro.1", "U_Astro.2", "U_Astro.3", "U_Astro.4") ~ "U_R_Astrocytes",
      receiver %in% c("U_Endo", "U_endo.ECs", "U_endo.FLCs", "U_endo.VSMC", "U_PER_END1", "U_PER_END2", "U_PER_END3") ~ "U_R_Endothelial",
      receiver %in% c("U_MG1", "U_MG2", "U_MG3", "U_Mic", "U_Micro", "U_Micro.0", "U_Micro.1", "U_Micro.2", "U_Micro.3", "U_Micro.4", "U_Micro.5", "U_Micro.6", "U_Micro.7", "U_Micro.8") ~ "U_R_Microglia",
      receiver %in% c("U_Ex1", "U_Ex0", "U_EX1", "U_EX2", "U_EX3", "U_EX4", "U_EX5", "U_Excit", "U_NeuronE.0", "U_NeuronE.1", "U_NeuronE.2", "U_NeuronE.6") ~ "U_R_NeuronE",
      receiver %in% c("U_In", "U_INH1", "U_INH2", "U_INH3", "U_INH4", "U_Inhit", "U_NeuronI.3", "U_NeuronI.4", "U_NeuronI.5") ~ "U_R_NeuronI",
      receiver %in% c("U_OPC", "U_OPC.0", "U_OPC.1", "U_OPC.2", "U_OPC.3", "U_OPC.4", "U_OPC.5", "U_OPC.6", "U_OPC1", "U_OPC2") ~ "U_R_OPC",
      receiver %in% c("U_ODC1", "U_ODC10", "U_ODC11", "U_ODC12", "U_ODC13", "U_ODC2", "U_ODC3", "U_ODC4", "U_ODC5", "U_ODC6", "U_ODC7", "U_ODC8", "U_ODC9",
                      "U_Oli0", "U_Oli1", "U_Oligo", "U_Oligo.0", "U_Oligo.1", "U_Oligo.2", "U_Oligo.3", "U_Oligo.4", "U_Oligo.5", "U_Oligo.6", "U_Oligo.7", "U_Oligo.8") ~ "U_R_Oligos",
      
      TRUE ~ receiver
    )
  ) %>% 
  select(-rank, -interaction, -mean_exp, -interaction_id, -interaction_id2) %>% 
  distinct()

gene_dictionary <- readRDS("~/crosstalk_webserver/data/genes_and_interactions.rds")
#gene_dictionary <- readRDS("data/genes_and_interactions.rds")



my_collapse <- function(x) {
  paste0(unique(x), collapse = ",")
}

# Interaction dictionary is the dataset where we have filtered the data from realdata to have ligand_pair and ligand and receptor details which are required for interaction table or gene table when search is performed
inter_dic <- realdata %>% 
  select(id_cp_interaction:is_integrin) %>% 
  select(-annotation_strategy) %>% 
  distinct() %>% 
  mutate(partner_a = gsub("simple:|complex:", "", partner_a)) %>% 
  mutate(partner_b = gsub("simple:|complex:", "", partner_b)) %>% 
  left_join(gene_dictionary, by = c("id_cp_interaction", "partner_a" = "key")) %>% 
  select(-hgnc_symbol, -ensembl, -partner, -starts_with("protein"), -uniprot, -source, -annotation_strategy) %>% 
  left_join(gene_dictionary, by = c("id_cp_interaction", "partner_b" = "key")) %>% 
  select(-hgnc_symbol, -ensembl, -partner, -starts_with("protein"), -uniprot) %>% 
  group_by(id_cp_interaction) %>% 
  summarise(across(everything(), my_collapse))
colnames(inter_dic) <- gsub("\\.x$", "_a", colnames(inter_dic))
colnames(inter_dic) <- gsub("\\.y$", "_b", colnames(inter_dic))



inter_dic <- inter_dic %>%
  mutate(
    receptor_a = as.logical(receptor_a),
    ligand = case_when(
      receptor_a & complex_name_b != "NA" ~ complex_name_b,
      receptor_a & complex_name_b == "NA" ~ gene_name_b,
      !receptor_a & complex_name_a != "NA" ~ complex_name_a,
      !receptor_a & complex_name_a == "NA" ~ gene_name_a,
      TRUE ~ NA_character_
    ),
    receptor = case_when(
      receptor_a & complex_name_a != "NA" ~ paste0(complex_name_a, ifelse(!is.na(gene_name_a), paste0(" {", gene_name_a, "}"), "")),
      receptor_a & complex_name_a == "NA" ~ gene_name_a,
      !receptor_a & complex_name_b != "NA" ~ paste0(complex_name_b, ifelse(!is.na(gene_name_b), paste0(" {", gene_name_b, "}"), "")),
      !receptor_a & complex_name_b == "NA" ~ gene_name_b,
      TRUE ~ NA_character_
    ),
    
    ligand_pair = paste(ligand, receptor, sep = " and ")
  )
inter_dic$ligand_pair <- gsub(" {.*?}", "", inter_dic$ligand_pair, perl = TRUE)






# Define the assigned colors for each categorySystematic analysis of cellular crosstalk reveals a role for SEMA6D-TREM2 regulating microglial function in Alzheimer’s disease

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#F781BF")
celltypes <- c("Astrocytes", "Microglia", "NeuronE", "NeuronI", "OPC", "Oligos", "Endothelial")

category_color <- setNames(colors, celltypes)