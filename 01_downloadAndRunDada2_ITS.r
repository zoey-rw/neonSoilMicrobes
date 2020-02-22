
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/trimPrimers.R")

library(dada2)
library(ShortRead)
library(Biostrings)
source("/projectnb/talbot-lab-data/zrwerbin/NEFI_microbe/paths.r")
source("dada2_fwd.reads.R")
source("trimPrimers.R")

library(foreach)
library(doParallel)
cores=detectCores()
registerDoParallel(cores)

amplicon <- "ITS"




# RUN DADA2 TO INFER EXACT SEQUENCE VARIANTS (ESVs)

in.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/filt_seqs/ITS/"
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/seq_tables/ITS"
seqRuns <- sort(list.files(in.dir, pattern = "....."))
tic()
for (s in 7:length(seqRuns)) {
  #for (s in 1:1) {
  tic()
  cat(paste("Running dada2 inference for sequencing run:", seqRuns[s]))
  filtpath <- file.path(in.dir, seqRuns[[s]])
  if (length(list.files(filtpath))==0) next()
  out.file <- paste0(out.dir, "/", seqRuns[s],".rds")
  runDada2_fwd.reads(filtpath, out.file)
  cat(paste("Exact sequence variant table constructed for ",  seqRuns[s],". Saving...\n"))
  toc()
}

# COMBINE ESV TABLES FROM EACH RUN, remove chimeras.
#load tables and merge, set output path.----
seq.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/seq_tables/ITS/"
seqRuns <- list.files(seq.dir, full.names = T)
#seqtab.nochim <- mergeRemoveChimeras(seqRuns = seqRuns, output.path = "data/output_files/ITS/mergedOtuTable.rds")
seqtab.nochim <- readRDS("data/output_files/ITS/mergedOtuTable.rds")

# ASSIGN TAXONOMY

tic()
unite_path     <- '/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data_old/sh_general_release_dynamic_01.12.2017.fasta'
to_assign <- colnames(seqtab.nochim)
tax.out <- dada2::assignTaxonomy(to_assign, unite_path, verbose = TRUE, multithread = T)
toc()
saveRDS(tax.out, "data/output_files/ITS/taxTable_legacy_fwdOnly.rds")
cat('Taxonomy output saved.\n')



abun_ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_groupAbundances_ITS.rds")


plot.list <- list()
for (group in c("Ectomycorrhizal","Arbuscular","Saprotroph","Pathogen")){
  samp.df <- sample_data(abun_ITS$ps_ITS)
  samp.df$asDate <- as.Date(samp.df$dates, "%Y%m%d")
  y <- abun_ITS[[group]]$rel.abundances[,group]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .75) + ggtitle(group) 
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])


