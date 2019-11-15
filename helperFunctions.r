# helper functions


# from the dada2 tutorial
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


# check for primers
checkPrimers <- function(fwd, rev, FWD.orients, REV.orients){
  print(
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fwd), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rev), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fwd), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = rev))
)
}

# run dada core inference algorithm, from dada2 big data tutorial
runDada2 <- function(filtpath, out.dir, ...){
  tic()
  # File parsing
  filtFs <- list.files(filtpath, pattern="_F_filt.fastq.gz", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="_R_filt.fastq.gz", full.names = TRUE)
  
  filtFs <- filtFs[sapply(strsplit(basename(filtFs), "_[RF]"), `[`, 1) %in% sapply(strsplit(basename(filtRs), "_[RF]"), `[`, 1) ]
  filtRs <- filtRs[sapply(strsplit(basename(filtRs), "_[RF]"), `[`, 1) %in% sapply(strsplit(basename(filtFs), "_[RF]"), `[`, 1) ]
  sample.names <- sapply(strsplit(basename(filtFs), "_[RF]"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  sample.namesR <- sapply(strsplit(basename(filtRs), "_[RF]"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  set.seed(100)
  # Learn forward error rates
  errF <- learnErrors(filtFs, nbases=1e6, multithread=TRUE)
  # Learn reverse error rates
  errR <- learnErrors(filtRs, nbases=1e6, multithread=TRUE)
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  for(i in 1:length(sample.names)) {
    sam <- sample.names[[i]]
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
    cat(paste("Exact sequence variants inferred for sample:",  sam,". \n"))
  }
  rm(derepF); rm(derepR)
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, out.dir)
  cat(paste("Saving sequence table at:", out.dir,". \n"))
  toc()
}



# tic/toc functions from Colin Averill:
#' Two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
tic <- function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc <- function() print(Sys.time()-timer)




#normalize otu table. function from Colin Averill:
pro.function <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}

#' Set truncation length
#' 
#' Decides on truncation length for trimmed reads based on quality score means. Default cutoff is a score of 30. 
#' Warns you if the beginning (first 10) bases are low-quality, but returns the first low-quality base after the 10th base.
#' If no bases are below score, returns the last base.
#' 
#'  fl : input files (prints suggested length for each; if used in script, could take the first input, or the lowest)
#'  qscore : default = 30.
#'  n : default = 5e+05 don't know exactly why this matters, but it's part of how 'qa' is called within the dada2 scripts...
#' 
#' 
get_truncation_length <-function (fl, qscore = 30, n = 5e+05, quiet = TRUE){
  trunc_lengths <- data.frame(file = character(0), early_lowqual = numeric(0),
                              trunc_length = numeric(0))
  for (f in fl) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    lowqual <- which(means < qscore)
    lowqual_in_first_20 <- length(which(lowqual <= 20))
    lowqual_after_20 <- lowqual[lowqual > 20][1]
    
    trunc_lengths <- rbind(trunc_lengths, data.frame(file = f, early_lowqual = lowqual_in_first_20,
                                                     trunc_length = ifelse(is.na(lowqual_after_20), length(means), lowqual_after_20))) 
    
    if (quiet == FALSE) {
      if(lowqual_in_first_20 > 0){
        cat(paste0(basename(f),': ', lowqual_in_first_20, ' bases within first 20 bases are below your quality score. Consider trimming left.\n'))
      } 
      if (is.na(lowqual_after_20)){
        cat(paste0(basename(f), ': After first 20 bases, no bases with a mean under your quality score found. Truncate at end of read, base: ', length(means),'\n'))
        
      } else if (!is.na(lowqual_after_20)){
        cat(paste0(basename(f),': After first 20 bases, first mean below your quality score is base: ',lowqual_after_20,'\n'))
        #return(lowqual_after_20)
        
      } else "Something's wrong. Inspect function/reads."
    } # end printing to console  
  } # end loop
  return(trunc_lengths$trunc_length)  
} # end function
