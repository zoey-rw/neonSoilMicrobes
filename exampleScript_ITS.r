source("helperFunctions.r")

source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("dada2_fwd.reads.R")
source("trimPrimers.R")

library(foreach)
library(doParallel)
cores=detectCores()
registerDoParallel(cores)








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





ps_ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_ITS.rds")
fg_ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/fg_abun.rds")
abun_ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_groupAbundances_ITS.rds")


plot.list <- list()
for (group in c("Ectomycorrhizal","Arbuscular","Saprotroph","Pathogen")){
  samp.df <- sample_data(abun_ITS$phylum$phyloseq)
  samp.df$asDate <- as.Date(samp.df$dates, "%Y%m%d")
 # samp.df <- samp.df[!grepl("2017", samp.df$dateID),]
  y <- fg_ITS$rel.abundances[,group]
 # y <- fg_ITS$rel.abundances[,group][which(!grepl("2017", samp.df$dateID))]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID)) + ggtitle(group)
}


