# start of function to run dada2 inference pipeline on NEON soil microbial samples
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/helperFunctions.r")

library(dada2)
library(ShortRead)
library(Biostrings)

in.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/filtSeqs/"
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/seqTab/"

if (!require(dada2)) stop("Please install dada2 according to the directions at: https://benjjneb.github.io/dada2/dada-installation.html")

seqRuns <- list.files(in.dir)

if (!dir.exists(out.dir)) dir.create(out.dir)

runDada2 <- function(filtpath, out.dir, ...){
  tic()
  # File parsing
  filtFs <- list.files(filtpath, pattern="_F_filt.fastq.gz", full.names = TRUE)
  filtRs <- list.files(filtpath, pattern="_R_filt.fastq.gz", full.names = TRUE)
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
  for(sam in sample.names) {
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
  saveRDS(seqtab, out.file)
toc()
}



tic()
for (s in 1:length(seqRuns)) {
  cat(paste("Running dada2 inference for sequencing run:", seqRuns[[s]]))
  filtpath <- file.path(in.dir, seqRuns[[s]])
  out.file <- paste0(out.dir, "/", seqRuns[s],".rds")
  runDada2(filtpath, out.file)
  cat(paste("Exact sequence variant table constructed for ",  seqRuns[s],". Saving...\n"))
}
toc()