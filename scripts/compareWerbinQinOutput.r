# phys <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10086.soil_phys.rds")
# phys <- phys
library(phyloseq)
library(dplyr)
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")

abun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_groupAbundances_16S.rds")

werbin_ps_orig <- abun$ps_16S
werbin_samp.df <- sample_data(werbin_ps_orig)

# phys$soilSampleID <- phys$sampleID
# phys$sampleID <- NULL
# df <- as(sample_data(abun$ps_16S), "data.frame")
# rownames(df) <- sample_names(abun$ps_16S)
# df$soilSampleID <- df$sample
# df$sampleID <- NULL
# df$dateID <- NULL
# df$row.names  <- rownames(df)
# phys <- phys[phys$soilSampleID %in% df$soilSampleID,]
# 
# samp.df <- merge(df, phys, all.x = T, all.y = F, by = "soilSampleID")
# samp.df <- samp.df[!duplicated(samp.df$row.names),]
# rownames(samp.df) <- samp.df$row.names
# samp.df$sampleID <- samp.df$row.names
# samp.df$legacy <- ifelse(samp.df$asDate.x < "2015-01-01", T, F)
# sample_data(ps) <- samp.df


qin_ps_orig <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/qin_NEON_16S_phyloseq_DL08-13-2019.Rds")
qin_samp.df <- sample_data(qin_ps_orig)
qin_samp.df$dateID <- substr(qin_samp.df$collectDate.x, 1, 7)

samp.keep <- intersect(werbin_samp.df$sample, qin_samp.df$sampleID)
qin_ps <- prune_samples(qin_samp.df$sampleID %in% samp.keep, qin_ps_orig)
werbin_ps <- prune_samples(werbin_samp.df$sample %in% samp.keep, werbin_ps_orig)

sample_names(qin_ps) <- make.unique(as.character(qin_samp.df$sampleID))
colnames(tax_table(qin_ps)) <- tolower(rank_names(qin_ps))


qin_phy <- get_tax_level_abun(qin_ps,tax_rank_list = c("phylum"), min_seq_depth = 50)
qin_phy_rel <- qin_phy$phylum$rel.abundances
colnames(qin_phy_rel) <- tolower(colnames(qin_phy_rel))
qin_phy_rel$sampleID <- rownames(qin_phy_rel)

zoey_phy_rel <- abun$phylum$rel.abundances
colnames(zoey_phy_rel) <- paste0(colnames(zoey_phy_rel), "_zoey")
zoey_phy_rel$sampleID <- gsub("-GEN-DNA.", "", rownames(zoey_phy_rel))

merged <- merge(zoey_phy_rel, qin_phy_rel, by = "sampleID")

plot(merged$acidobacteria_zoey ~ merged$acidobacteria)
plot(merged$actinobacteria_zoey ~ merged$actinobacteria)
