# takes the taxonomy and ESV tables from dada2, creates a phyloseq object, adds functional groups,
# and calculates the abundances/prevalence at every level.
library(phyloseq)
library(data.table)

# set output paths
otu_tax_output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/master_otu_tax.rds"
ps.output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_ITS.rds"
abun.output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_groupAbundances_ITS.rds"

# load recent data
otu_ITS_recent <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/mergedOtuTable.rds")
tax_ITS_recent <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/taxTable_ITS.rds")
# load legacy data
otu_ITS_legacy <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/otuTable_legacy_fwdOnly.rds")
tax_ITS_legacy <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/taxTable_legacy_fwdOnly.rds")

# remove extra sites from legacy data.
otu_ITS_legacy <- otu_ITS_legacy[which(substr(rownames(otu_ITS_legacy),1,4) %in% c("CPER","DSNY","HARV","OSBS","STER")),]
otu_ITS_legacy <- otu_ITS_legacy[,which(colSums(otu_ITS_legacy)>0)]
otu <- dada2::mergeSequenceTables(otu_ITS_recent, otu_ITS_legacy) # combine OTU tables

# combine taxonomy tables and reformat
tax <- rbind(tax_ITS_recent, tax_ITS_legacy)
for (i in 1:ncol(tax)) tax[, i] <- substring(tax[, i], 4) # remove leading "k__"
colnames(tax) <- tolower(colnames(tax))
tax <- as.data.frame(apply(tax,2,tolower)) # make everything lower case

# remove taxa that do not assign to a kingdom, and samples with less than 1k reads.
tax <- tax[tax$kingdom == 'fungi',] 
otu <- otu[rowSums(otu) >= 1000,]
tax <- tax[intersect(rownames(tax),colnames(otu)),]
otu <- otu[,intersect(rownames(tax),colnames(otu))]
identical(rownames(tax),colnames(otu))

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

#### ADD FUNCTIONAL GROUP ASSIGNMENTS ####
#Assign functions to taxa using FUNGuild. takes just a minute.
n <- detectCores()
n <- 8
tax.fun <- fg_assign_parallel(tax,n.cores=n)

#Build phylo-functional group taxonomy tables.----
fg <- data.table(tax.fun)
groups <- c('Arbuscular','Animal Pathogen','Plant Pathogen','Saprotroph','Wood Saprotroph','Ectomycorrhizal')
fg[,groups] <- 0
for (g in 1:length(groups)) {
  group <- groups[[g]]
  fg <- fg %>% 
    dplyr::mutate(!!group := dplyr::case_when(grepl(!!group, guild) ~ !!group,
                                              TRUE ~ "other"))
}
fg <- fg[,c('kingdom','phylum','class','order','family','genus','species', groups)]

# replace white spaces with underscores.
fg <- apply(fg, 2, function(x) gsub('\\s+', '_',x))
colnames(fg) <- lapply(colnames(fg), function(x) gsub('\\s+', '_',x))

# reassign to phyloseq object.
rownames(fg) <- taxa_names(ps)
tax_table(ps) <- fg

# get abundances for taxonomic and functional groups.
abun <- get_tax_level_abun(ps, tax_rank_list = rank_names(ps)[c(2:6,8:13)])
out <- c(ps_ITS, abun)
names(out) <- c("ps_ITS", names(abun))

# save output!!
saveRDS(ps, ps.output.path)
saveRDS(out, abun.output.path)
