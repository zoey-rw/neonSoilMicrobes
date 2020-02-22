# source functions, load libraries
library(dplyr)
source("helperFunctions.r")
source("dada2_fwd.reads.R")
source("trimPrimers.R")
source("qualityFilter.r")
source("mergeSequenceTables.r")
# set parameters for analysis
amplicon <- "16S"

#### CREATE OUTPUT LOCATIONS
# Output for raw sequences
raw.seqs <- "data/raw_seqs/16S/raw_from_neon/fastq/"
seqRuns <- sort(list.files(raw.seqs))
# set output directory for filtered sequences (separated by sequencing run)
filt.seqs <- "data/raw_seqs/16S/raw_from_neon/filt_seqs/"
# set output directory for sequence tables (separated by sequencing run)
seq.tabs <- "data/raw_seqs/16S/raw_from_neon/seq_tables/"

#### QUALITY-FILTER ####
track.file <- file.path(filt.seqs, paste0("trackFilt_", amplicon, ".txt")) # tracking file
out <- list()
for(s in 1:length(seqRuns)) {
  raw.seq.dir <- file.path(raw.seqs, seqRuns[s])  # trimmed file directory for this run
  out.dir <- file.path(filt.seqs, seqRuns[s])  # quality-filtered file directory for this run
  cat(paste("\n\nQuality-filtering sequencing run:", seqRuns[s],"\n\n\n"))
  
  reads <- list.files(raw.seq.dir, full.names = T)
  sample.names <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(reads))
  filt <- file.path(filt.seqs, seqRuns[s], "/", paste0(sample.names, "_filt.fastq.gz"), fsep ='')
  
  out[[s]] <- filterAndTrim(reads, filt, truncLen=250, #truncation.length,
                                   maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                                   compress=TRUE, multithread=TRUE)
  names(out)[[s]] <- seqRuns[[s]]
}
all.out <- rbind.named.dfs(out)
all.out <-  all.out %>% 
  mutate(filtered = (reads.out/reads.in)) %>%  
  group_by(seqRun) %>% 
  mutate(avg_filtered_seqRun = mean(filtered))
#rownames(all.out) <- rownames(rbind.named.dfs(out))
all.out$sampleID <- rownames(rbind.named.dfs(out))
write.table(all.out, track.file)



# RUN DADA2 TO INFER EXACT SEQUENCE VARIANTS (ESVs)
# Currently, dada2 is not using pooling (or pseudo-pooling) just because of how this code is set up,
# But changing this would improve detection of rare taxa 
seqRuns <- list.files(filt.seqs, "run......") #exclude tracking files
for (s in 1:length(seqRuns)) {
  cat(paste("Running dada2 inference for sequencing run:", seqRuns[[s]]))
  filtpath <- file.path(filt.seqs, seqRuns[[s]])
  out.file <- paste0(seq.tabs, "/", seqRuns[s],".rds")
  runDada2_fwd.reads(filtpath, out.file, seed = 100)
  cat(paste("Exact sequence variant table constructed for ",  seqRuns[[s]],". Saving...\n"))
}

# MERGE OTU TABLES FROM EACH RUN 
seqRuns <- list.files(seq.tabs, full.names = T)
seqtab.nochim <- mergeRemoveChimeras(seqRuns = seqRuns, output.path = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/raw_from_neon/seq_tables/legacyOtuTable_16S.rds")

# ASSIGN TAXONOMY  
greengenes.path <- "/projectnb/talbot-lab-data/NEFI_data/gg_13_8_train_set_97.fa"
tax_output_path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/raw_from_neon/seq_tables/legacyTaxTable_16S.rds"

seqtab.nochim <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/raw_from_neon/seq_tables/legacyOtuTable_16S.rds")
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

