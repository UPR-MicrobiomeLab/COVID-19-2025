#build phyloseq object from QIIME2 data with manual taxonomy curation
#comments in english only
#LOAD
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(OTUtable)
library(reshape2)
library(ape)
library(qiime2R)
library(tidyverse)
library(readr)
library(tidyr)
library(PERFect)
library(devtools)
library(remotes)
library(QsRutils)

setwd("/Users/danielavargasrobles/Library/CloudStorage/GoogleDrive-danielavargasrobles@gmail.com/My Drive/Filipa_Daniela_2021 -present/Proyectos/ORAL_COVID")

tree <- read_tree("QIIME2/rooted_tree.nwk")

phy <- qza_to_phyloseq(features = "QIIME2/table_filtered_wo5.qza")
ASV <- otu_table(phy)
data<-as.data.frame(t(ASV))


taxonomy <- read_tsv('QIIME2/taxonomy.tsv', show_col_types = FALSE)

taxonomy_cleaned0 <- taxonomy %>%
  dplyr::select(-Confidence) %>% # Eliminar la columna Confidence
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ gsub("^\\s+", "", .))) %>%
  mutate(across(everything(), ~ gsub(" ", "_", .))) %>%  # Reemplazar espacios por guiones bajos
  mutate(
    Kingdom = ifelse(is.na(Kingdom) | Kingdom == '' | Kingdom == 'k__', 'k__Unclassified', Kingdom),
    Phylum = ifelse(is.na(Phylum) | Phylum == '' | Phylum == 'p__', 'p__Unclassified', Phylum),
    Class = ifelse(is.na(Class) | Class == '' | Class == 'c__', 'c__Unclassified', Class),
    Order = ifelse(is.na(Order) | Order == '' | Order == 'o__', 'o__Unclassified', Order),
    Family = ifelse(is.na(Family) | Family == '' | Family == 'f__', 'f__Unclassified', Family),
    Genus = ifelse(is.na(Genus) | Genus == '' | Genus == 'g__', 'g__Unclassified', Genus),
    Species = ifelse(is.na(Species) | Species == '' | Species == 's__', 's__Unclassified', Species),
    Species = ifelse(Genus != 'g__Unclassified' & Species == 's__Unclassified', paste0(Genus, "_s__Unclassified"), 
                     ifelse(Genus != 'g__Unclassified' & Species != 's__Unclassified', paste0(Genus, " ", sub('s__', '', Species)), 
                            "g__Unclassified_s__Unclassified"))
  ) %>%
  column_to_rownames(var = "Feature ID") %>%  # Establecer Feature ID como nombres de fila
  filter(Kingdom != "k__Archaea" & Kingdom != "Unassigned" & Phylum != "p__Unclassified", Class != "c__Unclassified", Order != "o__Unclassified")


write_csv(taxonomy_cleaned0, 'QIIME2/taxonomy_cleaned.csv', quote = "none")

manual_taxonomy <- data.frame(
  ASV_ID = c(
    "45b6647181e654120901aa7ce1cdff2c", "d1ff3b3af7c615dd2d8355ec92f3ab94", 
    "d7ffc279960f1c0b98e8eb0cae19143e", "700414b6f32ae1482dd1b75ab4789759", 
    "64709e2311a97bcacdb19c400182c519", "e8f63d2c25f9f03db40ef17d0cafa54e", 
    "c2d1b01fb113a31b9aedd821a48dabd1", "90acb141c9db37919de356dd87320402", 
    "5065ee431ecd4d40b67abc23af2710b5", "a0e599dd22189b50ab3ba7dcb38ad388", 
    "3e93cbd3d9b47bf0c96c7bcbfc49b99e", "669b52db3057999401c70166473910eb", 
    "497c4ce5804e46f74e506bc48cefd6c2", "a20c343808b9a59137a3ebe7b53a59a3", 
    "4af9096c96276a92acf64e719fb066ce", "93ce1df4279f5632946eaa3acb062188", 
    "343fcdd806fecc1010a13ae5bb61f782", "121041944368eb7b79d9ea0e9b6d288c", 
    "d9200f9d884faec9320e36d626d9c067", "6885cbab9ea7d2c1478ab053af984914", 
    "f954084122b893907636934399e518d7", "a46071b81db4ced3c1ec35bf92b0c01c", 
    "6f3da9a32524f8e93ab4d71f480a3e00", "34d2624becb07cbce6f1a4e6f816651e", 
    "b9acce8396584e035a7b54d2e688bf0e", "cb57004a30385ade6e84b9176e7973cc", 
    "3405abb98df1d7f215329f13781c21e6", "af1a367211c27028b4d00cabc25a3df4", 
    "1c8175400f30ce0368598d5fefb30ff4"
  ),
  Genus = c(
    "g__Streptococcus", "g__Haemophilus", "g__Streptococcus", "g__Gemella", 
    "g__Neisseria", "g__Granulicatella", "g__Streptococcus", "g__Fusobacterium", 
    "g__Prevotella", "g__Prevotella", "g__Streptococcus", "g__Lancefieldella", 
    "g__Veillonella", "g__Neisseria", "g__Lautropia", "g__Oribacterium", 
    "g__Haemophilus", "g__Haemophilus", "g__Actinomyces", "g__Parvimonas", 
    "g__Streptobacillus", "g__Haemophilus", "g__Granulicatella", "g__Neisseria", 
    "g__Abiotrophia", "g__Parvimonas", "g__Mogibacterium", "g__Neisseria", 
    "g__Selenomonas"
  ),
  Species = c(
    "oralis", "parainfluenzae", "parasanguinis", "Unclassified", 
    "mucosa", "adiacens", "australis", "pseudoperiodonticum", 
    "melaninogenica", "jejuni", "DTU_2020_1000888_1_SI_GRL_NUU_041A", 
    "parvula", "parvula", "mucosa", "mirabilis", "sinus", 
    "influenzae", "parainfluenzae", "naeslundii", "Unclassified", 
    "moniliformis", "sputorum", "Unclassified", "flava", 
    "defectiva", "micra", "pumilum", "elongata", "Unclassified"
  )
)


# Add ASV_ID as a column in taxonomy_cleaned0
taxonomy_cleaned1 <- taxonomy_cleaned0 %>%
  rownames_to_column(var = "ASV_ID")

# Merge with manual taxonomy corrections
taxonomy_cleaned2 <- taxonomy_cleaned1 %>%
  left_join(manual_taxonomy, by = "ASV_ID", suffix = c("", ".manual")) %>%
  mutate(
    Genus = ifelse(!is.na(Genus.manual), Genus.manual, Genus),
    Species = ifelse(!is.na(Species.manual), Species.manual, Species)
  ) %>%
  dplyr::select(-Genus.manual, -Species.manual)

# Edit the Species column
taxonomy_cleaned3 <- taxonomy_cleaned2 %>%
  mutate(
    Species = ifelse(
      Species == "Unclassified",
      paste0(Genus, "_Unclassified"),  # Incluir el género antes de Unclassified
      ifelse(!grepl("^g__", Species), paste0(Genus, " ", Species), Species)
    ),
    Species = gsub(" ", "_", Species)  # Reemplaza espacios por guiones bajos en el nombre de la especie
  )

taxonomy_cleaned <- taxonomy_cleaned3 %>%
  column_to_rownames(var = "ASV_ID")

print(tail(taxonomy_cleaned))

# Convert taxonomy_cleaned to tax_table to phyloseq objects
taxonomy_cleaned <- tax_table(as.matrix(taxonomy_cleaned))

#filtered_otu_table <- otu_table(t(filtered_data_simultaneous), taxa_are_rows = TRUE)
filtered_otu_table <- otu_table(t(data), taxa_are_rows = TRUE)

# PRUNE taxonomy after ASV table filtering
#common_taxa <- intersect(taxa_names(taxonomy_cleaned), taxa_names(filtered_otu_table)) #descomentar si uso "perfect" for filter
common_taxa <- intersect(taxa_names(taxonomy_cleaned), taxa_names(data))
TAXA_clean <- prune_taxa(common_taxa, taxonomy_cleaned)
OTU_clean <- prune_taxa(common_taxa, filtered_otu_table)

# Create the phyloseq object
phy1 <- phyloseq(TAXA_clean, OTU_clean)
phy1


# 
# #####verify the unclassified at species level ASVs#####
# 
# # Filtrar las ASVs que no están clasificadas a nivel de especie
# unclassified_asvs <- tax_table(phy1) %>% 
#   as.data.frame() %>%
#   filter(grepl("Unclassified", Species)) %>%
#   rownames()
# 
# # Extraer la tabla de OTU del objeto phyloseq
# otu_table <- otu_table(phy1)
# 
# # Filtrar la tabla de OTU para quedarse solo con las ASVs sin clasificar en especie
# unclassified_otu_table <- otu_table[unclassified_asvs, ]
# 
# # Extraer la información taxonómica (familia y género)
# unclassified_taxonomy <- tax_table(phy1) %>%
#   as.data.frame() %>%
#   filter(rownames(.) %in% unclassified_asvs) %>%
#   dplyr::select(Family, Genus)
# 
# # Calcular el número total de lecturas (reads) para cada ASV sin clasificar
# unclassified_counts <- rowSums(unclassified_otu_table)
# 
# # Calcular el número de muestras que contienen cada ASV sin clasificar
# unclassified_sample_counts <- apply(unclassified_otu_table, 1, function(x) sum(x > 0))
# 
# # Crear un dataframe con los resultados, incluyendo la familia y el género
# unclassified_summary <- data.frame(
#   ASV_ID = unclassified_asvs,
#   Family = unclassified_taxonomy$Family,
#   Genus = unclassified_taxonomy$Genus,
#   Total_Reads = unclassified_counts,
#   Sample_Count = unclassified_sample_counts
# )
# 
# # Ordenar los resultados por el número total de lecturas (descendente)
# unclassified_summary %>% 
#    arrange(desc(Total_Reads))
# 















# Filtering steps

phy_filtered <- prune_taxa(taxa_sums(phy1) > 10, phy1) 

# Calculate the presence of each ASV across samples
asv_presence <- apply(otu_table(phy_filtered), 1, function(x) sum(x > 0))

# Filter ASVs that are present in more than 1 sample
phy1 <- prune_taxa(asv_presence > 1, phy_filtered)

# Extract the filtered OTU table
filtered_otu_table <- otu_table(otu_table(phy1), taxa_are_rows = TRUE)


common_taxa_filtered <- taxa_names(phy1) 

# Prune the taxonomy table to keep only the filtered ASVs
TAXA_clean_filtered <- prune_taxa(common_taxa_filtered, TAXA_clean)

# Generate the final phyloseq object with filtered data
PHY_filtered <- phyloseq(
  tax_table(as.matrix(TAXA_clean_filtered)), 
  otu_table(as.matrix(otu_table(phy1)), taxa_are_rows = TRUE), 
  phy_tree(tree)
)



#### Add metadata ####
met <- read.csv("metadata/met_13oct22_10_clean_25ago25_final.csv", header = TRUE, row.names = 1)


#Impute age
# First, calculate the mean age excluding the 0 values
media_edad <- met %>%
  dplyr::filter(Age > 0) %>% # Excluir las edades que son 0
  summarise(mean_age = mean(Age, na.rm = TRUE)) %>% # Calcular la media
  pull(mean_age) # Extraer el valor de la media

# Now, replace the 0 values in the Age column with the calculated mean age
met <- met %>%
  mutate(Age = ifelse(Age == 0, media_edad, Age))

# Select only relevant columns
met %>% 
  dplyr::select(Age, BMI, Sex, antibiotics_last_2months, Patient_ID, Pregnant)

# Make a new study_group3 variable
met=met %>%
  mutate(study_group3 = case_when(
    Pregnant == "Yes"& Sex=="Female" ~ "Pregnant",        # Si la persona está embarazada
    Sex == "Female"& Pregnant == "No" ~ "Female",            # Si el sexo es femenino
    Sex == "Male" & Pregnant == "No" ~ "Male"                 # Si el sexo es masculino
  ))                                   
table(met$study_group3)

# Reorder metadata to match the samples in the filtered phyloseq object
met_ordered_filtered <- met[match(sample_names(PHY_filtered), met$id), ]

# assign row names to the metadata
rownames(met_ordered_filtered) <- sample_names(PHY_filtered)

# Convert the ordered metadata to sample_data format
met_ordered_sample_data <- sample_data(met_ordered_filtered)

# Generate the final phyloseq object with metadata
PHYraw <- phyloseq(
  tax_table(tax_table(PHY_filtered)), 
  otu_table(otu_table(PHY_filtered), taxa_are_rows = TRUE), 
  sample_data(met_ordered_sample_data), 
  phy_tree(tree)
)


save(PHYraw,file="phy_objects/PHYraw.RData")
#load("phy_objects/PHYraw.RData")

#### Rarefaction #####

# Calcular la suma de secuencias por muestra
sumatoria <- sample_sums(PHYraw)
min(sumatoria)

sumatoria_ordenada <- sort(sumatoria)
print(head(sumatoria_ordenada, 30))

# Agregar la suma de secuencias como una nueva columna en la metadata
sample_data(PHYraw)$num_secuencias <- sumatoria

# Mostrar la metadata con la nueva columna
print(sample_data(PHYraw))

seq_to_check=data.frame(sample_data(PHYraw))
class(seq_to_check)
po=seq_to_check %>% 
  dplyr::select(num_secuencias, covid_casos_control_analysis, COVID_consensus) %>% 
  arrange(num_secuencias)
head(po, 40)
#write.csv(seq_to_check, "metadata/met_13oct22_9_clean_24april24_casos_control_column_nseqs.csv", row.names = TRUE)


#Rarefaction#######
# este valor de 3450, 9 muestras covid negativas todas menos de 28 sequencias
set.seed(711)
PHYrar=rarefy_even_depth(PHYraw, sample.size = 3418,#2142,#5430,#3489,
                         rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


#One note that can be useful, I also use the estimation of Good's indexes (library(QsRutils)) to evaluate, for instance, if the rarefaction approach is adequate, by looking at the median and mean values
#you can see that both mean and median values dropped a bit after rarefaction. If they remain >90% its ok otherwise the rarefied samples are not reflecting adequately the original taxa diversity
not_rare=summary(goods(otu_table(PHYraw)))
not_rare
rare=summary(goods(otu_table(PHYrar)))
rare

#save rarefied object
save(PHYrar,file="phy_objects/PHYrar.RData")
met=sample_data(PHYrar)

##ALPHA diversity indices#####
#1. extract ASV table 
OTU1 = as(otu_table(PHYrar), "matrix")
OTUdf = as.data.frame(OTU1)

#Calculating 4 different alpha diversity metrics
#at asv level
shannon = vegan::diversity(OTUdf, index = "shannon", MARGIN = 2, base = exp(1))
simpson = vegan::diversity(OTUdf, index = "simpson", MARGIN = 2)
#chao = apply(OTUdf, 2, chao1)
chao=estimate_richness(PHYrar, split = TRUE,measures = "Chao1")
obs = apply(OTUdf, 2, function(x){
  nws = length(which(x>0))
  return(nws)
})

div_asv = data.frame(shannon_rar_asv = shannon, 
                     simpson_rar_asv = simpson,
                     chao_rar_asv = chao,
                     obs_rar_asv = obs)



#merging columns of alpha metrics with the rest of the metadata by SampleID
md1 = sample_data(PHYrar)
md0=as.data.frame(as(md1, "matrix"))
write.csv(md0, "md0.temp", row.names = TRUE)
md = read.csv("md0.temp")

asv=cbind(md,div_asv)
rownames(asv) = sample_names(PHYrar)

sampledata = sample_data(asv)

#Re doing phyloseq object with new added columns (alpha div metrics) in metadata
PHYrar_alpha = phyloseq(otu_table(PHYrar),tax_table(PHYrar), sampledata, phy_tree(PHYrar) )
hola=sample_data(PHYrar_alpha)
hola$shannon_rar_asv
save(PHYrar_alpha,file="phy_objects/PHYrar_alpha.RData")
#load("phy_objects/PHYrar_alpha.RData")

