
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
library(dada2)
library(ShortRead)
library(Biostrings)

in.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/ITS/trimmedSeqs"
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/ITS/filtSeqs"
log.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/logfiles/ITS"

seqRuns <- list.files(in.dir)

# read in raw sequences
tic()
out.list <- list()
for(s in 1:length(seqRuns)){
  
  trimmed.run.dir <- file.path(in.dir, seqRuns[s])  # read in primer-trimmed file directory for each run
  filtered.run.dir <- file.path(out.dir, seqRuns[s]) # create filtered file directory for each run
  dir.create(filtered.run.dir, suppressWarnings(paste0("In dir.create(filtered.run.dir) : '", 
                                                       filtered.run.dir, "' already exists")))
  fnFs <- sort(list.files(trimmed.run.dir, pattern = "_R1.fastq.gz", full.names = TRUE)) # primer-trimmed FWD files
  fnRs <- sort(list.files(trimmed.run.dir, pattern = "_R2.fastq.gz", full.names = TRUE)) # primer-trimmed REV files
  filtFs <- file.path(filtered.run.dir, basename(fnFs)) # filtered FWD files
  filtRs <- file.path(filtered.run.dir, basename(fnRs)) # filtered REV files
  
  if(length(filtFs) != length(filtRs)) {
    cat(paste0("Forward and reverse files do not match for: ",seqRuns[s])); next()
  }
  if(length(fnFs)==0 | length(fnRs)==0) {
    cat(paste0("Missing all forward or all reverse reads for: ",seqRuns[s])); next()
  }
  
  #### FILTERING AND TRIMMING ENDS OF READS ####
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 50,
                       maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE, matchIDs = TRUE,
                       compress=TRUE, verbose=TRUE, multithread=12)
  out <- cbind(seqRuns[[s]], "quality_filter", out)
  
  # saving each tracking file for primer trimming
  write.table(out, file.path(log.dir, paste0("filtTrack",seqRuns[[s]],".txt")))
  out.list[[s]] <- out
}

# combining into one tracking file to save
filtTrack.all <- do.call(rbind, out.list)
write.table(filtTrack.all, file.path(log.dir, "filtTrack_allITSruns.txt"))
toc()
