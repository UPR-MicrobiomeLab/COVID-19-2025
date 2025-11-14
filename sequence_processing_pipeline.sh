# University of Puerto Rico
# N files: 1216
# 16S V4 sequences
# QIIME2 DADA2

EXPLANATION  
The analysis consists of:  
1. Processing the sequences using Qiime2, which is described in this file  
2. The data analysis, done in R using a phyloseq object. This object was generated from the ASV table, the rooted tree, the taxonomy table (from Qiime), and the metadata.  
The generation of the phyloseq object from these 4 files is in "principal_oral_vag_cov.R", and the phyloseq object already generated is in "/phy_objects/"

# I downloaded the files from "box" using Firefox.  
# I transferred them to my hard drive and uncompressed them there, and from there I ran the Qiime2 import into ".qza"  

# activate QIIME2
conda activate qiime2-2021.11

# remove all files except
find . \! -name '*.gz' -delete

# count files
find *.gz -type f | wc -l

# File COVID_Filipa_2022_2 on the external hard drive

# I downloaded the files and placed them on the hard drive; from there I ran the manifest.
# I used the manifest method because the other method didn’t allow me to upload it.
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/COVID_19_filipa/README/manifest2.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# It was run again, removing these samples from the manifest because they were duplicates
# (Filipa later confirmed this based on what a student told her on WhatsApp)
COV-T1-616  
COV-T1-616  
COV-T3-040  
COV-T3-927  
COV-T3-040  
COV-T3-927

# This sample appeared this way in the metadata, but I changed it in the metadata to COV-T1-584
# so it would match its file name.
COV-T1-590-584

# Sequence quality summary
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization summary_demux.qzv

qiime tools view summary_demux.qzv

# DADA2 algorithm
qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-demux.qza \
--p-trim-left-f 10 \
--p-trunc-len-f 249 \
--p-trim-left-r 10 \
--p-trunc-len-r 239 \
--o-representative-sequences dada2_paired-end-demux.qza \
--o-denoising-stats stat \
--o-table table.qza

# Visualization of ASV table
qiime metadata tabulate \
  --m-input-file table.qza \
  --o-visualization table.qzv
qiime tools view table.qzv

#################
# DATA TRAINING #
#################

# Importing database
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/VMC_gg_homd_14_5_db.fa \
  --output-path /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/gg_ext.qza

# Importing taxonomy table associated with database
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/VMC_gg_homd_14_5_taxonomy_no_duplicatesDV_formated.txt \
  --output-path /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/gg_ext_tax.qza

# Extract reference reads from the database; for this we use V4 primers
qiime feature-classifier extract-reads \
  --i-sequences /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/gg_ext.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-min-length 240 \
  --p-max-length 255 \
  --o-reads /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/ref-seqs_gg_ext_V4_240_255bp.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/ref-seqs_gg_ext_V4_240_255bp.qza \
  --i-reference-taxonomy /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/gg_ext_tax.qza \
  --o-classifier /Users/danielavargasrobles/Google\ Drive/Filipa/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/classifier.qza

#####################################
# Now we run it against the data    #
#####################################

# classified sequences
qiime feature-classifier classify-sklearn \
  --i-classifier  /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/VAGINAL/DATABASE_FOR_VAGINAL_MB/CVM_db_16S/classifier.qza \
  --i-reads dada2_paired-end-demux.qza \
  --o-classification taxonomy.qza

# FROM TABLE: Removing Mitochondria and Chloroplast
# This script is case-insensitive, so "mitochondria" returns the same as "Mitochondria".
# Only taxa classified up to phylum were kept, since "k__Bacteria" was not informative.
qiime taxa filter-table \
--i-table table.qza \
--i-taxonomy taxonomy.qza  \
--p-include p__ \
--p-exclude Mitochondria,Chloroplast,cyanobacteria \
--o-filtered-table table_filtered.qza

# Visualization of filtered ASV table
qiime metadata tabulate \
  --m-input-file table_filtered.qza \
  --o-visualization table_filtered.qzv
qiime tools view table_filtered.qzv

# These samples are in the ASV table but not in the metadata;
# however, I already included them in met_13oct_2.csv
COV-NEGATIVE  
COV-T3-511  
COV-T3-539  
COV-T2-532  
COV-T1-590

qiime feature-table filter-samples \
  --i-table table_filtered.qza \
  --m-metadata-file /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/COVID_19_filipa/Metadata/met_13oct22_8.tsv \
  --p-where "[id] IN ('COV-NEGATIVE','COV-T3-511','COV-T3-539','COV-T2-532','COV-T1-590')" \
  --p-exclude-ids TRUE \
  --o-filtered-table table_filtered_wo5.qza

# Visualization of ASV table
qiime metadata tabulate \
  --m-input-file table_filtered_wo5.qza \
  --o-visualization table_filtered_wo5.qzv
qiime tools view table_filtered_wo5.qza.qzv

# barplot
# To load the metadata I tried many methods to generate a .tsv file but failed.
# I had to convert it using this CSV→TSV converter:
# https://products.groupdocs.app/conversion/csv-to-tsv
qiime taxa barplot \
  --i-table table_filtered_wo5.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/COVID_19_filipa/Metadata/met_13oct22_8.tsv \
  --o-visualization taxa_barplot.qzv

qiime tools view taxa_barplot.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences dada2_paired-end-demux.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table_filtered_wo5.qza \
  --p-sampling-depth 4000 \
  --m-metadata-file /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/COVID_19_filipa/Metadata/met_13oct22_8.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/COVID_19_filipa/Metadata/met_13oct22_8.tsv \
  --o-visualization core-metrics-results/shannon_vector-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file /Volumes/GoogleDrive/My\ Drive/Filipa_Daniela_2021/COVID_19_filipa/Metadata/met_13oct22_8.tsv \
  --m-metadata-column COVID-Status \
  --o-visualization core-metrics-results/bray-body-site-significance.qzv \
  --p-pairwise

qiime tools view core-metrics-results/bray_curtis_emperor.qzv  
qiime tools view core-metrics-results/bray-body-site-significance.qzv  
qiime tools view core-metrics-results/shannon_vector-group-significance.qzv  
qiime tools view core-metrics-results/bray_curtis_emperor.qzv

qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported-tree

qiime tools export \
  --input-path taxonomy.qza

