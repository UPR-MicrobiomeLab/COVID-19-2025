# COVID-19-2025
## SARS-CoV-2 Oral Microbiome Analysis

This repository contains the code, phyloseq objects, metadata, and QIIME 2 pipeline used to process and analyze the oral microbiome data for the manuscript. All microbiome analyses were performed in R using phyloseq, and sequence processing was performed in QIIME 2.

### Repository Structure

`sequence_processing_pipeline.sh` — QIIME 2 processing pipeline.

`scripts/1_generate_phyloseq_obj.R` — Builds the main phyloseq object.
`scripts/2_analyses_and_plots.Rmd` — Runs all analyses and generates the plots included in the paper.

`data/phy_object/PHYrar_alpha.RData` — Main phyloseq object (all available samples).

`data/phy_object/PHYrar_alpha_to_paper.rds` — Filtered phyloseq object (only the samples included in the paper).

`data/raw files/taxonomy.tsv` — Taxonomy table exported from QIIME 2.

`data/raw files/table_filtered_wo5.qza` — Feature table (QIIME 2).

`data/raw files/rooted_tree.nwk` — Phylogenetic tree.

`data/raw files/met_13oct22_10_clean_25ago25_final.csv` — Raw metadata file used to generate the phyloseq object.

`data/metadata/metadata_PHYrar_alpha_to_paper.csv` — Metadata containing only the samples included in the final manuscript.

## Phyloseq Objects and Usage

Two phyloseq objects are provided under `data/phy_object/:

`PHYrar_alpha.RData` (Main phyloseq object used by all scripts): Includes all samples available after preprocessing. This is the object that is loaded and used throughout the scripts by default.

`PHYrar_alpha_to_paper.rds` (A filtered version of the main object): Includes only the samples that were actually used in the manuscript. Use this object if you want your analyses to match exactly the sample set in the paper.

Recommendation:

To run or explore all scripts as written → use `PHYrar_alpha.RData`

To use only the samples set used in the paper → use `PHYrar_alpha_to_paper.rds`

## Metadata

The metadata exclusivelly with the samples included in the manuscript is stored in: `data/metadata/metadata_PHYrar_alpha_to_paper.csv`. It also includes three helper columns with Yes / No values to indicate which samples were used for each of the three main analyses:

`sarscov2_pos_vs_neg_crosssectional`
→ Cross-sectional evaluation of the influence of SARS-CoV-2 infection on oral microbiome.

`sarscov2_pos_vs_neg_longitudinal`
→ Longitudinal evaluation of the influence of SARS-CoV-2 infection on the oral microbiome.

`vaccination_influence_longitudinal`
→ Longitudinal evaluation of the influence of vaccination on the oral microbiome.

These columns can be used directly in R to filter the phyloseq object or metadata to the subset of samples relevant to each analysis.

## Scripts
`scripts/1_generate_phyloseq_obj.R`: This script generates the main phyloseq object `PHYrar_alpha.RData` using the raw inputs in data/raw files/:

`table_filtered_wo5.qza`

`taxonomy.tsv`

`rooted_tree.nwk`

`met_13oct22_10_clean_25ago25_final.csv`


The script imports, cleans, and merges these inputs to create the phyloseq object stored in `data/phy_object/PHYrar_alpha.RData`.

`scripts/2_analyses_and_plots.Rmd`: This R notebook document runs the analyses and creates the figures used in the manuscript, including:

1. Alpha diversity analyses

2. Beta diversity and ordinations

3. PERMANOVA and related statistical models

4. Plots and summary outputs

By default, it loads the main phyloseq object: `load("data/phy_object/PHYrar_alpha.RData")`

If you want to restrict analyses strictly to the samples included in the paper, you can instead load: `phy <- readRDS("data/phy_object/PHYrar_alpha_to_paper.rds")` and adapt the R notebook accordingly.

##QIIME 2 Processing Pipeline

sequence_processing_pipeline.sh: contains the QIIME 2 command-line pipeline used to process the raw sequencing data, including:

1. Import of raw reads 

2. Denoising and ASV generation

3. Taxonomy assignment

4. Filtering steps
 
5. Export of feature table, taxonomy, and tree (`table_filtered_wo5.qza`, `taxonomy.tsv`, `rooted_tree.nwk`)

These outputs are then used as inputs to build the phyloseq object in `1_generate_phyloseq_obj.R`.
