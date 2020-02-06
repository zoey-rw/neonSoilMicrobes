
runDada2 <- function(filtpath, out.file, seed = NULL, ...){
  tic()
  if (!dir.exists(out.dir)) dir.create(out.dir)
  if (!is.null(seed)) set.seed(seed)

  # File parsing
  filtFs <- list.files(filtpath, pattern="_R1.fastq.gz", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="_R2.fastq.gz", full.names = TRUE)
  
  filtFs <- filtFs[sapply(strsplit(basename(filtFs), "_R"), `[`, 1) %in% sapply(strsplit(basename(filtRs), "_R"), `[`, 1) ]
  filtRs <- filtRs[sapply(strsplit(basename(filtRs), "_R"), `[`, 1) %in% sapply(strsplit(basename(filtFs), "_R"), `[`, 1) ]
  sample.names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  sample.namesR <- sapply(strsplit(basename(filtRs), "_R"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
  
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  # Learn forward error rates
  errF <- learnErrors(filtFs, nbases=1e7, multithread=TRUE)
  # Learn reverse error rates
  errR <- learnErrors(filtRs, nbases=1e7, multithread=TRUE)
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
    merger <- mergePairs(ddF, derepF, ddR, derepR, maxMismatch=1)
    mergers[[i]] <- merger
    cat(paste("Exact sequence variants inferred for sample:",  sam,". \n"))
  }
  #rm(derepF); rm(derepR)
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, out.file)
  toc()
}