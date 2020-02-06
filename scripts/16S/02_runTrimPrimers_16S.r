# this script trims primers and quality-filters all non-legacy 16S data from NEON.
rm(list=ls())

source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/trimPrimers.R")

library(dada2)
library(ShortRead)
library(Biostrings)

# read in raw sequences
raw.seqs <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/fastq/"
# set output directory for trimmed sequences (separated by sequencing run)
all.out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/trimmed_seqs/16S/"
seqRuns <- sort(list.files(raw.seqs))
amplicon <- "16S"
tic()

 for(s in 1:length(seqRuns)) {
  cat(paste("\n\nTrimming primers and filtering reads for sequencing run:", seqRuns[s],"\n\n\n"))

#### CREATING FILE NAMES/DIRECTORIES ####
raw.seq.dir <- file.path(raw.seqs, seqRuns[s])  # raw file directory
out.dir <- file.path(all.out.dir, seqRuns[s])  # primer-trimmed file directory
track.file <- file.path(all.out.dir, "log", paste0("trackPrimers_", seqRuns[s], "_", amplicon, ".csv"))

  ### TRIMMING PRIMERS ####
    trimPrimers(raw.seq.dir = raw.seq.dir,
                amplicon = "16S",
                out.dir = out.dir,
                FWD <- "CCTACGGGNBGCASCAG",  #341F primer, listed in PCR amplification file
                REV <- "GACTACNVGGGTATCTAATCC", #Pro805R, this rev is listed in the PCR amplification file
                flip_rev_primer = F, remove.filtN = T,
                test = F, track.file = track.file,
                trimLeft = 0)
}

