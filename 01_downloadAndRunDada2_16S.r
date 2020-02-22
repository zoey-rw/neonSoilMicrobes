# source functions, load libraries
library(dplyr)
source("helperFunctions.r")
source("dada2.R")
source("trimPrimers.R")
source("qualityFilter.r")
source("mergeSequenceTables.r")
# set parameters for analysis
amplicon <- "16S"
wd <- getwd() # change this path if you want data to download elsewhere

#### CREATE OUTPUT LOCATIONS
# Creates output folders for each step of analysis
create_data_directory(path = wd) # creates a "data" directory with subdirectories
# Output for raw sequences
raw.seqs <- "data/raw_seqs/16S/fastq/"
seqRuns <- sort(list.files(raw.seqs))
# set output directory for trimmed sequences (separated by sequencing run)
trimmed.seqs <- "data/trimmed_seqs/16S/"
# set output directory for filtered sequences (separated by sequencing run)
filt.seqs <- "data/filt_seqs/16S/"
# set output directory for sequence tables (separated by sequencing run)
seq.tabs <- "data/seq_tables/16S/"

#### DOWNLOAD RAW SEQUENCES ####
samples <- downloadRawSequenceData(site = "HARV", amplicon = amplicon, 
                        startdate="2017-05", enddate="2017-06", 
                        out.dir = "/data/raw_seqs/ITS", check.size=F)

#### TRIM PRIMERS ####
track.file <- file.path(trimmed.seqs, paste0("trackPrimers_", amplicon, ".txt")) # tracking file
out <- list()
for(s in 1:length(seqRuns)) {
  raw.seq.dir <- file.path(raw.seqs, seqRuns[s])  # raw file directory for this run
  out.dir <- file.path(trimmed.seqs, seqRuns[s])  # primer-trimmed file directory for this run
  cat(paste("\n\nTrimming primers for sequencing run:", seqRuns[s],"\n\n\n"))
  out[[s]] <- trimPrimers(raw.seq.dir = raw.seq.dir,
                          #samples = samples,
                          out.dir = out.dir)
  names(out)[[s]] <- seqRuns[[s]]
}
all.out <- rbind.named.dfs(out)
write.table(all.out, track.file)


#### QUALITY-FILTER ####
track.file <- file.path(filt.seqs, paste0("trackFilt_", amplicon, ".txt")) # tracking file
out <- list()
for(s in 1:length(seqRuns)) {
  #for(s in 1:1) {
  if(seqRuns[[s]]=="BRPH4") next()
    trimmed.seq.dir <- file.path(trimmed.seqs, seqRuns[s])  # trimmed file directory for this run
  out.dir <- file.path(filt.seqs, seqRuns[s])  # quality-filtered file directory for this run
  cat(paste("\n\nQuality-filtering sequencing run:", seqRuns[s],"\n\n\n"))
  out[[s]] <- qualityFilter(trimmed.seq.dir = trimmed.seq.dir,
                          #samples = samples,
                          out.dir = out.dir)
  names(out)[[s]] <- seqRuns[[s]]
}
  all.out <- rbind.named.dfs(out)
  all.out <-  all.out %>% 
    mutate(filtered = (reads.out/reads.in)) %>%  
    group_by(seqRun) %>% 
    mutate(avg_filtered_seqRun = mean(filtered))
  rownames(all.out) <- rownames(rbind.named.dfs(out))
  write.table(all.out, track.file)
  
  
  
# RUN DADA2 TO INFER EXACT SEQUENCE VARIANTS (ESVs)
  # Currently, dada2 is not using pooling (or pseudo-pooling) just because of how this code is set up,
  # But changing this would improve detection of rare taxa 
  seqRuns <- list.files(filt.seqs, "\\b.....\\b") #exclude tracking files
for (s in 1:length(seqRuns)) {
  cat(paste("Running dada2 inference for sequencing run:", seqRuns[[s]]))
  filtpath <- file.path(filt.seqs, seqRuns[[s]])
  out.file <- paste0(seq.tabs, "/", seqRuns[s],".rds")
  runDada2(filtpath, out.file, seed = 100)
  cat(paste("Exact sequence variant table constructed for ",  seqRuns[[s]],". Saving...\n"))
}

# MERGE OTU TABLES FROM EACH RUN 
seqRuns <- list.files(seq.tabs, full.names = T)
seqtab.nochim <- mergeRemoveChimeras(seqRuns = seqRuns, output.path = "data/output_files/16S/otuTable_16S.rds")

# ASSIGN TAXONOMY  
greengenes.path <- "/projectnb/talbot-lab-data/NEFI_data/gg_13_8_train_set_97.fa"
tax_output_path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/taxTable_16S.rds"



seqtab.nochim <- readRDS("data/output_files/16S/otuTable_16S.rds")
all.out <- read.table(track.file)
track <- cbind.data.frame(all.out, nonchim = rowSums(seqtab.nochim))
track <- track %>% 
           mutate(final_retained = nonchim/reads.in)

to_assign <- colnames(seqtab.nochim)
#assign taxonomy.
tic()
cat('Assigning taxonomy using the greengenes Classifier...\n')
out <- dada2::assignTaxonomy(to_assign, greengenes.path, multithread = T, verbose = TRUE)
cat('Taxonomy assignment complete. ')
toc()
saveRDS(out, tax_output_path)


# combine OTU and TAX into phyloseq object in make_ps_16S.r
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/createTaxFunction.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/addBacterialFunction.r")

# read in phyloseq, add function
ps_16S <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_16S.rds")
tax_fun_ref <-  createTaxFunction()
tax_df = as.data.frame(as(tax_table(ps_16S), "matrix"))
tax_with_function <- addBacterialFunction(tax = tax_df, tax_fun_ref = tax_fun_ref)
rownames(tax_with_function) <- taxa_names(ps_16S)
tax_table(ps_16S) <- as.matrix(tax_with_function)

abun <- get_tax_level_abun(ps_16S, tax_rank_list = rank_names(ps_16S)[c(2:6,8:20)])
out <- c(ps_16S, abun)
names(out) <- c("ps_16S",names(abun))

saveRDS(abun, "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/groupAbundances.rds")

library(phyloseq)
abun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/groupAbundances.rds")

rel_abun <- sapply(abun, "[[", 2)
rel_abun$df <- sample_data(abun$phylum$phyloseq)
saveRDS(rel_abun, "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/relAbundances.rds")
rel_abun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/relAbundances.rds")
rel_abun$df <- rel_abun$df[!grepl("2017", rownames(rel_abun$df)),]

p <- list()
for (group in names(rel_abun)[6:18]){
  rel_abun[[group]] <- rel_abun[[group]][rownames(rel_abun[[group]]) %in% rownames(rel_abun$df),]
  
  df <- cbind.data.frame(asDate = rel_abun$df$asDate,
                         siteID = rel_abun$df$siteID,
                         horizon = rel_abun$df$horizon,
                         y = rel_abun[[group]][,group])
  p[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = 0.5) + ggtitle(group)
}
gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]])
gridExtra::grid.arrange(p[[7]], p[[8]], p[[9]], p[[10]], p[[11]], p[[12]])

