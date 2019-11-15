# start of function to run dada2 inference pipeline on NEON soil microbial samples
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
tic()
library(dada2)
library(ShortRead)
library(Biostrings)
library(foreach)
library(doParallel)

# detect and register cores for parallel computing
cores=detectCores()
registerDoParallel(cores)

in.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/filtSeqs/"
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/seqTab/"


if (!require(dada2)) stop("Please install dada2 according to the directions at: https://benjjneb.github.io/dada2/dada-installation.html")
out.list <- list()
seqRuns <- list.files(in.dir)
if (!dir.exists(out.dir)) dir.create(out.dir)


#out.list <- foreach(s=1:length(seqRuns), .errorhandling='pass') %dopar% {
#  out.list <- foreach(s=11:19, .errorhandling='pass') %dopar% {
    result <- tryCatch({
      cat(paste("Running dada2 inference for sequencing run:", seqRuns[[s]]))
      filtpath <- file.path(in.dir, seqRuns[[s]])
      out.file <- paste0(out.dir, "/", seqRuns[s],".rds")
      runDada2(filtpath, out.file)
      #cat(paste("Exact sequence variant table constructed for ",  seqRuns[s],". Saving...\n"))
      x <- "test"
      return(x)
    }, error=function(e) paste0("ERROR inferring sequence table for: ", seqRuns[s],"; skipping to next sequencing run.")
    )
    
    return(result)
  }
  toc()
  