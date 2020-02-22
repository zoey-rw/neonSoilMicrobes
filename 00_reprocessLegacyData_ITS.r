# re-process legacy ITS data
library(dada2)
library(ShortRead)
library(Biostrings)

source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/mergeSequenceTables.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEFI_microbe/paths.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/dada2_fwd.reads.R")

library(foreach)
library(doParallel)
cores=detectCores()
registerDoParallel(cores)

amplicon <- "ITS"
seqRuns <- list()
seqRuns[[1]] <- NEON_ITS_run150922_r1_fastq.dir
seqRuns[[2]] <- NEON_ITS_run150225_r1_fastq.dir

#### CREATING FILE NAMES/DIRECTORIES ####
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/filt_seqs/ITS/150922"  # primer-trimmed file directory
filt.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/filt_seqs/ITS"
log.dir <- file.path(filt.dir, "log")

# No need to trim primers, primers were trimmed when Colin received sequences.

# quality filter:
# first legacy run
track.file <- file.path(log.dir, paste0("trackFilt_150922_", amplicon, ".csv"))
filtered.run.dir <- file.path(filt.dir, "150922") # create filtered file directory for each run
fnFs <- sort(list.files(NEON_ITS_run150922_r1_fastq.dir, pattern = ".fastq", full.names = TRUE)) # primer-trimmed FWD files
filtFs <- file.path(filtered.run.dir, basename(fnFs)) # filtered FWD files

#### FILTERING AND TRIMMING ENDS OF READS ####
out <- filterAndTrim(fnFs, filtFs, minLen = 50,
                     maxN=0, maxEE=5, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=T)
write.table(out, track.file)

# second legacy run
track.file <- file.path(log.dir, paste0("trackFilt_150225_", amplicon, ".csv"))
filtered.run.dir <- file.path(filt.dir, "150225") # create filtered file directory for each run
fnFs <- sort(list.files(NEON_ITS_run150225_r1_fastq.dir, pattern = ".fastq", full.names = TRUE)) # primer-trimmed FWD files
filtFs <- file.path(filtered.run.dir, basename(fnFs)) # filtered FWD files
#### FILTERING AND TRIMMING ENDS OF READS ####
out <- filterAndTrim(fnFs, filtFs, minLen = 50,
                     maxN=0, maxEE=5, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=T)
write.table(out, track.file)

# run dada2
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/seq_tables/ITS"
  tic()
  cat(paste("Running dada2 inference for sequencing run: 150922"))
  filtpath <- file.path(filt.dir, "150922") 
  out.file <- paste0(out.dir, "/150922.rds")
  runDada2_fwd.reads(filtpath, out.file)
  cat(paste("Exact sequence variant table constructed for 150922. Saving...\n"))
toc()

  tic()
  cat(paste("Running dada2 inference for sequencing run: 150225"))
  filtpath <- file.path(filt.dir, "150225") 
  out.file <- paste0(out.dir, "/150225.rds")
  runDada2_fwd.reads(filtpath, out.file)
  cat(paste("Exact sequence variant table constructed for 150225. Saving...\n"))
  toc()

  
  # MERGE OTU TABLES FROM EACH RUN 
  seqtab.nochim <- mergeRemoveChimeras(seqRuns = c(paste0(out.dir, "/150922.rds"),
                                                   paste0(out.dir, "/150225.rds")), output.path = "data/output_files/ITS/otuTable_legacy_fwdOnly.rds")
  
  
  seqtab.nochim <- readRDS("data/output_files/ITS/otuTable_legacy_fwdOnly.rds")
  unite_path     <- '/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/reference_data/sh_general_release_dynamic_01.12.2017.fasta'
  
  tic()
  to_assign <- colnames(seqtab.nochim)
  tax.out <- dada2::assignTaxonomy(to_assign, unite_path, verbose = TRUE, multithread = T)
  toc()
  cat('Taxonomy assignment complete! yeahhhh.\n')
  saveRDS(tax.out, "data/output_files/ITS/taxTable_legacy_fwdOnly.rds")
  cat('Taxonomy output saved.\n')
  
  
