
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/trimPrimers.R")

library(dada2)
library(ShortRead)
library(Biostrings)
library(foreach)
library(doParallel)
cores=detectCores()
cores = 3
registerDoParallel(cores)

raw.seqs <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/ITS/fastq/"
all.out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/trimmed_seqs/ITS/"
seqRuns <- list.files(raw.seqs)
amplicon <- "ITS"
tic()
for (s in 1:length(seqRuns)) {
#  for (s in 3:6) {
  
    # skip runs I've already processed/checked
    if (s < 6) next()  
  
    cat(paste("\nTrimming primers and filtering reads for sequencing run:", seqRuns[s],"\n\n\n"))
    
    #### CREATING FILE NAMES/DIRECTORIES ####
    raw.seq.dir <- file.path(raw.seqs, seqRuns[s])  # raw file directory
    out.dir <- file.path(all.out.dir, seqRuns[s])  # primer-trimmed file directory
    track.file <- file.path(all.out.dir, "log", paste0("trackPrimers_", seqRuns[s], "_", amplicon, ".csv"))
    
    # Change trimLeft parameter (for dada2::filterAndTrim) because of low quality beginnings of reads
     if (seqRuns[s] %in% c("BTV4F", "C24VW", "C25T6","BRPH4")){
       cat("Trimming left 15 basepairs (run is either BTV4F, C24VW, or C25T6)")
       trimLeft <- 15
     } else trimLeft <- 0
   
    ### TRIMMING PRIMERS ####
    trimPrimers(raw.seq.dir = raw.seq.dir,
                amplicon = "ITS",
                out.dir = out.dir,
                FWD = "CTTGGTCATTTAGAGGAAGTAA",  #341F primer, listed in PCR amplification file
                REV = "GCTGCGTTCTTCATCGATGC", #Pro805R, this rev is listed in the PCR amplification file
                flip_rev_primer = F, remove.filtN = T,
                test = F, track.file = track.file,
                trimLeft = trimLeft)
}
cat("Primers trimmed for all runs.")
toc()