# not yet tested for ITS
# only uses forward reads for ITS

qualityFilter <- function(trimmed.seq.dir, out.dir, amplicon = "16S", fwdEE=5, revEE=5){
  fnFs <- sort(list.files(trimmed.seq.dir, pattern = "_R1.fastq.gz|_R1.fastq", full.names = TRUE)) # primer-trimmed FWD files
  fnRs <- sort(list.files(trimmed.seq.dir, pattern = "_R2.fastq.gz|_R2.fastq", full.names = TRUE)) # primer-trimmed REV files
  
  sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) # extract sample names
  filtFs <- file.path(out.dir, paste0(sample.names, "_filt_R1.fastq.gz")) # filtered FWD files
  filtRs <- file.path(out.dir, paste0(sample.names, "_filt_R2.fastq.gz")) # filtered REV files
  
  if(amplicon=="16S"){
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
  
  if (rev.trunc.length < 200) rev.trunc.length <- 200
  if (fwd.trunc.length < 245) fwd.trunc.length <- 245
  
  cat(paste0("Fwd truncation length: ", fwd.trunc.length, "\nRev truncation length: ", rev.trunc.length))
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(fwd.trunc.length, rev.trunc.length),
                       maxN=0, maxEE=c(fwdEE,revEE), truncQ=2, rm.phix=TRUE, matchIDs = TRUE,
                       compress=TRUE,
                       multithread=TRUE)
  
} else {
  #### FILTERING AND TRIMMING ENDS OF READS ####
  out <- filterAndTrim(fnFs, filtFs, minLen = 50,
                       maxN=0, maxEE=fwdEE, truncQ=2, rm.phix=TRUE,
                       compress=TRUE, verbose=TRUE, multithread=T)
}
  return(out)
}