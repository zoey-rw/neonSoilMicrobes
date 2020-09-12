# takes the taxonomy and ESV tables from dada2, creates a phyloseq object, adds functional groups,
# and calculates the abundances/prevalence at every level.
library(phyloseq)
library(data.table)
source("helperFunctions.r")
source("binTaxGroups.r")
source("addBacterialFunction.r")
source("createTaxFunction.r")
source("../NEFI_microbe/NEFI_functions/fg_assign_parallel.r")

# load recently-constructed data
otu_16S_recent <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/otuTable_16S.rds")
tax_16S_recent <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/taxTable_16S.rds")

# load legacy data
# otu_16S_legacy <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/NEON_processed/NEON_dada2_SV_table.rds")
# tax_16S_legacy <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/NEON_processed/NEON_dada2_tax_table.rds")

otu_16S_legacy <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/raw_from_neon/seq_tables/legacyOtuTable_16S.rds")
tax_16S_legacy <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/raw_from_neon/seq_tables/legacyTaxTable_16S.rds")

# set output paths
otu_tax_output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/master_otu_tax.rds"
ps.output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_16S.rds"
abun.output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_groupAbundances_16S.rds"

# Sample names are wrong for bacteria; this creates a table to link deprecatedVialID and geneticSampleID
map <- readRDS("/projectnb/talbot-lab-data/zrwerbin/pre-release-map_16S.rds")[,c("geneticSampleID","sampleID")]
map <- transform(map, row.names = sampleID, sampleID = NULL,
                 geneticSampleID = gsub('-GEN','',map$geneticSampleID))
map <- map[!duplicated(map),,drop=F]
rownames(otu_16S_legacy) <- gsub("_filt.fastq.gz", "", rownames(otu_16S_legacy))
otu_16S_legacy <- as.matrix(transform(merge(otu_16S_legacy,map,by=0), row.names=geneticSampleID, 
                                       Row.names=NULL, geneticSampleID=NULL))

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

abun <- get_tax_level_abun(ps, tax_rank_list = rank_names(ps)[c(2:6,8:20)])

out <- c(ps, abun)
names(out) <- c("ps_16S",names(abun))
saveRDS(out, abun.output.path)


# view functional groups
plot.list <- list()
for (group in names(abun)[7:19]){
  samp.df <- sample_data(abun$ps_16S)
  y <- abun[[group]]$rel.abundances[,group]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  #df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) 
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])
gridExtra::grid.arrange(plot.list[[5]], plot.list[[6]], plot.list[[7]], plot.list[[8]])
gridExtra::grid.arrange(plot.list[[9]], plot.list[[10]], plot.list[[11]], plot.list[[12]])
plot.list[[13]]



# view phyla
plot.list <- list()
prev <- abun$phylum$prevalence
phyla <- prev[order(-prev$prevalence),]
rel <- abun$phylum$rel.abundances
sort(apply(rel, 2, mean))
for (phylum in colnames(abun)[1:5]){
  samp.df <- sample_data(abun$ps_16S)
  y <- abun[[group]]$rel.abundances[,group]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  #df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) 
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])
gridExtra::grid.arrange(plot.list[[5]], plot.list[[6]], plot.list[[7]], plot.list[[8]])
gridExtra::grid.arrange(plot.list[[9]], plot.list[[10]], plot.list[[11]], plot.list[[12]])




# create a version with only >5000 reads

highreads <- prune_samples(rowSums(otu_table(abun$ps_16S)) > 6000, abun$ps_16S)


plot.list <- list()
for (group in names(abun)[7:19]){
  samp.df <- sample_data(highreads)
  y <- abun[[group]]$rel.abundances[,group]
  y <- y[which(rowSums(otu_table(abun$ps_16S)) > 6000)]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  #df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) 
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])
gridExtra::grid.arrange(plot.list[[5]], plot.list[[6]], plot.list[[7]], plot.list[[8]])
gridExtra::grid.arrange(plot.list[[9]], plot.list[[10]], plot.list[[11]], plot.list[[12]])
plot.list[[13]]





tax_df = as.data.frame(as(tax_table(abun_16S$ps_16S), "matrix"))
tax_fun_ref <-  createTaxFunction()
tax_fun_ref2 <- tax_fun_ref
tax_fun_ref2[tax_fun_ref2$taxon=="alphaproteobacteria",]$oligotroph <- 1
tax_with_function <- addBacterialFunction(tax = tax_df, tax_fun_ref = tax_fun_ref2)
rownames(tax_with_function) <- taxa_names(abun_16S$ps_16S)
tax_table(abun_16S$ps_16S) <- as.matrix(tax_with_function)

cop_olig <- get_tax_level_abun(abun_16S$ps_16S, tax_rank_list = rank_names(abun_16S$ps_16S)[c(19:20)])

cop_olig_plots <- list()
for (group in names(cop_olig)){
samp.df <- sample_data(abun_16S$ps_16S)
y <- cop_olig[[group]]$rel.abundances[,group]
df <- cbind.data.frame(asDate = samp.df$asDate,
                       siteID = samp.df$siteID,
                       horizon = samp.df$horizon,
                       y = y)
#df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
cop_olig_plots[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
  geom_point(aes(color = siteID)) +
  geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) + ylim(0, 1)
}

cop_olig_plots[[1]]
cop_olig_plots[[2]]
gridExtra::grid.arrange(cop_olig_plots[[1]], cop_olig_plots[[2]], nrow=2)

gridExtra::grid.arrange(plot.list[[13]], cop_olig_plots[[2]], nrow=1)
gridExtra::grid.arrange(plot.list[[12]]+ ylim(0, 1), cop_olig_plots[[1]], nrow=1)
