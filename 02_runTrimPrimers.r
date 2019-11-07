# this script trims primers and quality-filters all non-legacy 16S data from NEON.
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/trimPrimers.R")

library(dada2)
library(ShortRead)
library(Biostrings)

raw.seqs <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/fastq/"
parent.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/"
#raw.seqs <- list.files("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/fastq/B69PP/", full.names = T)
seqRuns <- list.files(raw.seqs)

# read in raw sequences
track <- list()
trimTrack <- list()
filtTrack <- list()
for (s in 15:length(seqRuns)) {
#  for (s in 2:2) {
    
  #### CREATING FILE NAMES/DIRECTORIES, PART 1 ####
  
  raw.run <- file.path(raw.seqs, seqRuns[s])  # raw file directory
  trimmed.run <- file.path(parent.dir, "trimmedSeqs", seqRuns[s])  # primer-trimmed file directory
  
  #### TRIMMING PRIMERS ####
  tryCatch({
  trim.qa <- trimPrimers(raw.seqs = raw.run, 
                         out.dir = trimmed.run, 
              FWD <- "CCTACGGGNBGCASCAG",  #341F primer, listed in PCR amplification file
              REV <- "GACTACNVGGGTATCTAATCC", #Pro805R, this rev is listed in the PCR amplification file
              flip_rev_primer = F, 
              test = FALSE)
  # 
  # raw.seqs = raw.run;
  # out.dir = trimmed.run; 
  # FWD <- "CCTACGGGNBGCASCAG";  #341F primer, listed in PCR amplification file
  # REV <- "GACTACNVGGGTATCTAATCC"; #Pro805R, this rev is listed in the PCR amplification file
  # flip_rev_primer = F;
  # test = FALSE;
  # remove.filtN = F;
  
  trimTrack[[s]] <- trim.qa
  names(trimTrack)[[s]] <- seqRuns[[s]]
  saveRDS(trimTrack[[s]], paste0("logfiles/trimTrack",seqRuns[[s]],".rds"))
  }, error=function(e){paste0("ERROR in trimming primers for: ", seqRuns[[s]],"; skipping to next sequencing run.")})

  #### CREATING FILE NAMES/DIRECTORIES, PART 2 ####
  
  fnFs <- sort(list.files(trimmed.run, pattern = "_R1.fastq", full.names = TRUE)) # primer-trimmed FWD files
  fnRs <- sort(list.files(trimmed.run, pattern = "_R2.fastq", full.names = TRUE)) # primer-trimmed REV files
  sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) # extract sample names
  filtered.run <- file.path(parent.dir, "filtSeqs", seqRuns[s]) # filtered file directory
  filtFs <- file.path(filtered.run, paste0(sample.names, "_F_filt.fastq.gz")) # filtered FWD files
  filtRs <- file.path(filtered.run, paste0(sample.names, "_R_filt.fastq.gz")) # filtered REV files
  
  #### GETTING TRUNCATION LENGTH USING QUALITY SCORES ####
  # just using the first 20. probably fine.
  
fwd.trunc.lengths <- list()
for (i in 1:10){
  if (file.exists(fnFs[i])){
  file.name <- fnFs[i]
  fwd.trunc.lengths[i] <- get_truncation_length(file.name, quiet=FALSE)
  } else next
}

rev.trunc.lengths <- list()
for (i in 1:10){
  if (file.exists(fnRs[i])){
  file.name <- fnRs[[i]]
  rev.trunc.lengths[i] <- get_truncation_length(file.name, quiet=FALSE)
  } else next
}

# remove any samples that have low-quality reads early on, to avoid spoiling the whole run
# idk how to loop this command
to_remove <- which(unlist(fwd.trunc.lengths) < 100 | unlist(rev.trunc.lengths) < 100)
if (length(to_remove)>0){
fwd.trunc.lengths <- fwd.trunc.lengths[-to_remove]
rev.trunc.lengths <- rev.trunc.lengths[-to_remove]
fnFs <- fnFs[-to_remove]
filtFs <- filtFs[-to_remove]
fnRs <- fnRs[-to_remove]
filtRs <- filtRs[-to_remove]
}

fwd.trunc.length <- round(min(unlist(fwd.trunc.lengths))) # tested this with mean vs minimum quality-score threshold;
rev.trunc.length <- round(min(unlist(rev.trunc.lengths))) # minimum retained more reads after filtering.

#### FILTERING AND TRIMMING ENDS OF READS ####

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(fwd.trunc.length, rev.trunc.length),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
colnames(out) <- c("beforeFilterReadCount", "afterFilterReadCount")
out <- as.data.frame(out)
out$run <- as.character(seqRuns[[s]])

# saving one tracking file for primer trimming
saveRDS(out, paste0("logfiles/trackTrim",seqRuns[[s]],".rds"))
}
# 
# # combining into one tracking file to save
# filtTrack.all <- do.call(rbind, filtTrack)
# saveRDS(filtTrack.all, "trackFilt.rds")
