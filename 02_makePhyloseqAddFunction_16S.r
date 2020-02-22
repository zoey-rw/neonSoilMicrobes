# takes the taxonomy and ESV tables from dada2, creates a phyloseq object, adds functional groups,
# and calculates the abundances/prevalence at every level.
library(phyloseq)
library(data.table)

# load recently-constructed data
otu_16S_recent <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/otuTable_16S.rds")
tax_16S_recent <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/taxTable_16S.rds")

# load legacy data
# otu_16S_legacy <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/NEON_processed/NEON_dada2_SV_table.rds")
# tax_16S_legacy <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/NEON_processed/NEON_dada2_tax_table.rds")

otu_16S_legacy <- t(readRDS("/projectnb/talbot-lab-data/caverill/to_retire/NEFI_microbial_data/map_otu/16S_otu_clean.rds"))
tax_16S_legacy <- readRDS("/projectnb/talbot-lab-data/caverill/to_retire/NEFI_microbial_data/map_otu/16S_tax_clean.rds")
colnames(tax_16S_legacy) <- c("kingdom","phylum","class","order","family","genus","species","ID")
tax_16S_legacy$ID <- NULL


# set output paths
otu_tax_output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/master_otu_tax.rds"
ps.output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_16S.rds"
abun.output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16/ps_groupAbundances_16S.rds"

# Sample names are wrong for bacteria; this creates a table to link deprecatedVialID and geneticSampleID
# map <- readRDS("/projectnb/talbot-lab-data/zrwerbin/pre-release-map_16S.rds")[,c("geneticSampleID","sampleID")]
# map <- transform(map, row.names = sampleID, sampleID = NULL, 
#                  geneticSampleID = gsub('-GEN','',map$geneticSampleID))
# map <- map[!duplicated(map),,drop=F]

# Remove extra sites from otu_ITS_legacy, combine OTU tables
otu_16S_legacy <- otu_16S_legacy[which(substr(rownames(otu_16S_legacy),1,4) %in% c("CPER","DSNY","HARV","OSBS","STER")),]
otu_16S_legacy <- otu_16S_legacy[,which(colSums(otu_16S_legacy)>0)]
rownames(otu_16S_legacy) <- make.unique(rownames(otu_16S_legacy))
otu <- dada2::mergeSequenceTables(otu_16S_recent, otu_16S_legacy)

# combine taxonomy tables and reformat
colnames(tax_16S_recent) <- tolower(colnames(tax_16S_recent))
tax <- rbind(tax_16S_recent, tax_16S_legacy)
for (i in 1:ncol(tax)) tax[, i] <- substring(tax[, i], 4) # remove leading "k__"
colnames(tax) <- tolower(colnames(tax))
tax <- as.data.frame(apply(tax,2,tolower)) # make everything lower case

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'bacteria',] # removes ~3000 ESVs
otu <- otu[rowSums(otu) >= 1000,]
tax <- tax[intersect(rownames(tax),colnames(otu)),]
otu <- otu[,intersect(rownames(tax),colnames(otu))]
identical(rownames(tax),colnames(otu))
rownames(otu) <- gsub("_16S_filt", "", rownames(otu))

#save this output.
saveRDS(list(otu,tax), otu_tax_output.path)

# create sample info and phyloseq.
df <- parseNEONsampleIDs(rownames(otu))
rownames(df) <- rownames(otu)

# put everything together into a phyloseq object
ps <- phyloseq(otu_table(otu, taxa_are_rows = FALSE), tax_table(as.matrix(tax)), sample_data(df))
# remove long DNA sequences - can be accessed with refseq(ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("SV", seq(ntaxa(ps)))

# save output!!
saveRDS(ps, ps.output.path)

#### ADD FUNCTIONAL GROUP ASSIGNMENTS ####
tax_fun_ref <-  createTaxFunction()
tax_df = as.data.frame(as(tax_table(ps), "matrix"))
tax_with_function <- addBacterialFunction(tax = tax_df, tax_fun_ref = tax_fun_ref)
rownames(tax_with_function) <- taxa_names(ps)
tax_table(ps) <- as.matrix(tax_with_function)

abun <- get_tax_level_abun(ps, tax_rank_list = rank_names(ps_16S)[c(2:6,8:20)])

out <- c(ps, abun)
names(out) <- c("otu",names(abun))
saveRDS(out, abun.output.path)





