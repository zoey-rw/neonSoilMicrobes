# this script trims primers and quality-filters all non-legacy 16S data from NEON.
rm(list=ls())

source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/trimPrimers.R")

library(dada2)
library(ShortRead)
library(Biostrings)
library(foreach)
library(doParallel)
cores=detectCores()
registerDoParallel(cores)

raw.seqs <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/raw_sequences_16S/recent/fastq/"
parent.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/"
seqRuns <- list.files(raw.seqs)

# # read in raw sequences
# out.list <- list()
# out.list <- foreach(s=1:length(seqRuns), .errorhandling='pass') %dopar% {
# #out.list <- list()
# #out.list <-  foreach(s=1:2) %dopar% {
#      #for (s in 1:2) {
#   cat(paste("\n\n\n\n\n\n\n\nTrimming primers and filtering reads for sequencing run:", seqRuns[s],"\n\n\n"))  
#   
#   #### CREATING FILE NAMES/DIRECTORIES, PART 1 ####
  # s <- 3
  # raw.run <- file.path(raw.seqs, seqRuns[s])  # raw file directory
  # trimmed.run <- file.path(parent.dir, "trimmedSeqs", seqRuns[s])  # primer-trimmed file directory
#   
#   ### TRIMMING PRIMERS ####
#   result <- tryCatch({
#   trimPrimers(raw.seq.dir = raw.run,
#                          out.dir = trimmed.run,
#               FWD <- "CCTACGGGNBGCASCAG",  #341F primer, listed in PCR amplification file
#               REV <- "GACTACNVGGGTATCTAATCC", #Pro805R, this rev is listed in the PCR amplification file
#               flip_rev_primer = F, remove.filtN = F,
#               test = FALSE, trackFile = TRUE)
# # 
  # raw.seq.dir = raw.run;
  # out.dir = trimmed.run;
  # FWD <- "CCTACGGGNBGCASCAG";  #341F primer, listed in PCR amplification file
  # REV <- "GACTACNVGGGTATCTAATCC"; #Pro805R, this rev is listed in the PCR amplification file
  # flip_rev_primer = F;
  # test = FALSE;
  # remove.filtN = F;
  # cutadapt.path = "/share/pkg.7/cutadapt/1.18/install/bin/cutadapt"
  # trackFile = TRUE
#  }, error=function(e){paste0("ERROR in trimming primers for: ", seqRuns[s],"; skipping to next sequencing run.")})
# 
#   return(result)
# }
# 




# read in raw sequences
out.list <- list()
out.list <- foreach(s=1:length(seqRuns), .errorhandling='pass') %dopar% {

  
result <- tryCatch({
  
    raw.run <- file.path(raw.seqs, seqRuns[s])  # raw file directory
    trimmed.run <- file.path(parent.dir, "trimmedSeqs", seqRuns[s])  # primer-trimmed file directory
  #### CREATING FILE NAMES/DIRECTORIES, PART 2 ####

  fnFs <- sort(list.files(trimmed.run, pattern = "_R1.fastq.gz", full.names = TRUE)) # primer-trimmed FWD files
  fnRs <- sort(list.files(trimmed.run, pattern = "_R2.fastq.gz", full.names = TRUE)) # primer-trimmed REV files

  sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) # extract sample names
  filtered.run <- file.path(parent.dir, "filtSeqs", seqRuns[s]) # filtered file directory
  filtFs <- file.path(filtered.run, paste0(sample.names, "_F_filt.fastq.gz")) # filtered FWD files
  filtRs <- file.path(filtered.run, paste0(sample.names, "_R_filt.fastq.gz")) # filtered REV files

  #### GETTING TRUNCATION LENGTH USING A SUBSET OF QUALITY SCORES ####
fwd.trunc.lengths <- list()
for (i in 1:10){
  if (file.exists(fnFs[i])) fwd.trunc.lengths[i] <- get_truncation_length(fnFs[i], #quiet=FALSE,
                                                                          qscore = 23)
   else next
}

rev.trunc.lengths <- list()
for (i in 1:10){
  if (file.exists(fnRs[i]))  rev.trunc.lengths[i] <- get_truncation_length(fnRs[[i]], #quiet=FALSE,
                                                                           qscore = 23)
   else next
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

if (seqRuns[[s]]=="BFDG8"){
fwd.trunc.length <- 250
rev.trunc.length <- 220
}
#### FILTERING AND TRIMMING ENDS OF READS ####

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(fwd.trunc.length, rev.trunc.length),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, matchIDs = TRUE,
                     compress=TRUE,
                     multithread=TRUE
)


# out <- filterAndTrim(fnFs[1:5], filtFs[1:5], fnRs[1:5], filtRs[1:5], truncLen=c(fwd.trunc.length, rev.trunc.length),
#                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, matchIDs = TRUE,
#                      compress=TRUE,
#                      multithread=TRUE
# )
colnames(out) <- c("beforeFilterReadCount", "afterFilterReadCount")
out <- as.data.frame(out)
out$run <- as.character(seqRuns[[s]])
# saving one tracking file for primer trimming
saveRDS(out, paste0("logfiles/trackFilt",seqRuns[[s]],".txt"))
out.list[[s]] <- out
return(head(out))
#
# # combining into one tracking file to save
}, error=function(e){paste0("ERROR in trimming primers for: ", seqRuns[s],"; skipping to next sequencing run.")}
  )

return(result)
}


filtTrack.all <- do.call(rbind, out.list)
saveRDS(filtTrack.all, "logfiles/trackFilt_all.rds")