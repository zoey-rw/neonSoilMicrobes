
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/trimPrimers.R")

library(dada2)
library(ShortRead)
library(Biostrings)
library(foreach)
library(doParallel)
cores=detectCores()
registerDoParallel(cores)

raw.seqs <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/raw_sequences_ITS/recent/fastq/"
parent.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/"
seqRuns <- list.files(raw.seqs)


for (s in 1:length(seqRuns)) {
    cat(paste("\n\n\n\n\n\n\n\nTrimming primers and filtering reads for sequencing run:", seqRuns[s],"\n\n\n"))
    
    #### CREATING FILE NAMES/DIRECTORIES, PART 1 ####
    raw.run <- file.path(raw.seqs, seqRuns[s])  # raw file directory
    trimmed.run <- file.path(parent.dir, "trimmedSeqs/ITS/", seqRuns[s])  # primer-trimmed file directory
    
    ### TRIMMING PRIMERS ####
    trimPrimers(raw.seq.dir = raw.run,
                out.dir = trimmed.run,
                FWD <- "CTTGGTCATTTAGAGGAAGTAA",  #341F primer, listed in PCR amplification file
                REV <- "GCTGCGTTCTTCATCGATGC", #Pro805R, this rev is listed in the PCR amplification file
                flip_rev_primer = F, remove.filtN = F,
                test = FALSE, trackFile = TRUE)
  }