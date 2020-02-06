
runDada2_fwd.reads <- function(filtpath, out.file, seed = NULL, ...){
  tic()
  if (!is.null(seed)) set.seed(seed)
  
  filtFs <- list.files(filtpath, pattern="fastq.gz|.fastq", full.names=TRUE) # CHANGE if different file extensions
  sample.names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1) # Assumes filename = samplename_RX.fastq.gz
  names(filtFs) <- sample.names
  # Learn error rates
  err <- learnErrors(filtFs, nbases = 1e7, multithread=TRUE, randomize=TRUE)
  # Infer sequence variants
  dds <- vector("list", length(sample.names))
  names(dds) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derep <- derepFastq(filtFs[[sam]])
    dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
    cat(paste("Exact sequence variants inferred for sample:",  sam,". \n"))
  }
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(dds)
  cat(paste0("\nDimensions of ESV table:\nSamples: ", dim(seqtab)[1], "\n\nESVs: ",dim(seqtab)[2]))
  saveRDS(seqtab, out.file)
  toc()
}
