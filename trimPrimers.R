# trim primers using cutadapt and the dada2 tutorial instructions.
# 
# # Input path to sequences:
# raw.seqs <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/fastq/"
# out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/trimmed/"
# # primer sequences come from dat$mmg_soilPcrAmplification file
# FWD <- "CCTACGGGNBGCASCAG"  
# REV <- "GTGYCAGCMGCCGCGGTAA" # don't know how i learned this, but this the 515F FWD primer for EMP.
# flip_rev_primer <- T
# REV <- "GACTACNVGGGTATCTAATCC" #this is the reverse primer shared by NEON, but the above one is actually found...
# flip_rev_primer <- F
# test <- TRUE
# cutadapt.path <- "/share/pkg.7/cutadapt/1.18/install/bin/cutadapt"  #input path to cutadapt

# example:
# trimPrimers(raw.seqs = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/fastq/", out.dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/trimmed/", flip_rev_primer = T, test = FALSE)

trimPrimers <- function(raw.seq.dir, out.dir,  samples = NULL, FWD = "CCTACGGGNBGCASCAG", REV = "GACTACNVGGGTATCTAATCC", cutadapt.path = "/share/pkg.7/cutadapt/1.18/install/bin/cutadapt", test = FALSE, remove.filtN = T, amplicon = "16S", trimLeft = 0, quiet=T){
  
  source("helperFunctions.r")
  tic()
  # Read in raw sequences
  fnFs <- sort(list.files(raw.seq.dir, pattern = "_R1.fastq.gz|_R1.fastq", full.names = TRUE, recursive = T))
  fnRs <- sort(list.files(raw.seq.dir, pattern = "_R2.fastq.gz|_R2.fastq", full.names = TRUE, recursive = T))
  
  if(!is.null(samples)) {
    fnFs <- fnFs[samples]
    fnRs <- fnRs[samples]
  }
  
  # Remove samples without forward and reverse reads 
  fnFs.names <- sapply(strsplit(basename(fnFs), split="_R"), '[', 1)
  fnRs.names <- sapply(strsplit(basename(fnRs), split="_R"), '[', 1)
  overlap <- intersect(fnFs.names, fnRs.names)
  fnFs <- fnFs[fnFs.names %in% overlap]
  fnRs <- fnRs[fnRs.names %in% overlap]
  
  if (length(overlap) == 0) {
    return("No overlap between forward and reverse reads.")
  }
  
  if (test){ # subset if testing
    fnFs <- head(fnFs, 1)
    fnRs <- head(fnRs, 1)
  }
  
  # Create output directory
  if(!dir.exists(out.dir)) dir.create(out.dir)
  
  # Create primer-trimmed file paths in output directory
  fnFs.cut <- file.path(out.dir, basename(fnFs))
  fnRs.cut <- file.path(out.dir, basename(fnRs))
  
  if(amplicon == "16S") {
    # Just trim the length of fwd and rev primers
    qa.out <- filterAndTrim(fnFs, fnFs.cut, fnRs, fnRs.cut, multithread = TRUE,
                              matchIDs = TRUE, trimLeft = c(nchar(FWD), nchar(REV)), compress=F)
   
     # Much more complicated for ITS...
    } else if (amplicon == "ITS") {
  
  # Create  N-filtered file paths in filtN/ subdirectory
  fnFs.filtN <- file.path(out.dir, "filtN", basename(fnFs)) 
  fnRs.filtN <- file.path(out.dir, "filtN", basename(fnRs))
  
  # First remove any N reads (Cutadapt/dada2 can't handle)
  out <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, truncQ = 0,
                       multithread = TRUE,
                       matchIDs = TRUE, trimLeft = trimLeft, compress=F)
  # Check for CutAdapt
  if (!is.null(cutadapt.path)) system(paste0("export PATH=",cutadapt.path,":$PATH"))
  if (system("which cutadapt") == 1) stop("Cannot find path to cutadapt. Set using the 'cutadapt.path' argument, or run from another location.")
  
  # Create reverse-complements of primers to search for
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  
  if(!quiet) {
    report <- "--quiet"
  } else report <- "--report=minimal"
  
  # Run Cutadapt
  for(i in 1:length(fnFs)) {
    system(paste("cutadapt ", report, R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads,
                 "-m 50", # minimum length of 50 reads 
                 "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                 fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
  
  # Remote N-filtered directory
  if(remove.filtN){
    cat("Removing temporary directory of N-filtered reads...")
    system(paste0("rm -r ", file.path(out.dir, "filtN")))
  }
  
  #### Create table to track reads through N-filter and primer trimming ####
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  if (length(fnFs.cut) < 2) { # just in case there's only one read.
    qaSummary.after.trim <- ShortRead::qa(fnFs.cut[[1]])[["readCounts"]][,1, drop=F]
    primers.before <- checkPrimers.wide(fnFs[[1]], fnRs[[1]], FWD.orients, REV.orients)
    names(primers.before) <- paste0("BeforeTrim_", names(primers.before))
    primers.after <- checkPrimers.wide(fnFs.cut[[1]], fnRs.cut[[1]], FWD.orients, REV.orients)
    names(primers.after) <- paste0("AfterTrim_", names(primers.after))
    qa.out <- cbind(out[1,, drop=F], qaSummary.after.trim, primers.before, primers.after)
  } else {
  # Get counts after trimming primers
  qaSummary.after.trim <- ShortRead::qa(fnFs.cut[1:2])[["readCounts"]][,1, drop=F]
  # Get counts for primer hits before trimming
  primers.before <- rbind(checkPrimers.wide(fnFs[[1]], fnRs[[1]], FWD.orients, REV.orients),
                          checkPrimers.wide(fnFs[[2]], fnRs[[2]], FWD.orients, REV.orients))
  names(primers.before) <- paste0("BeforeTrim_", names(primers.before))
  # Get counts for primer hits after trimming
  primers.after <- rbind(checkPrimers.wide(fnFs.cut[[1]], fnRs.cut[[1]], FWD.orients, REV.orients),
                         checkPrimers.wide(fnFs.cut[[2]], fnRs.cut[[2]], FWD.orients, REV.orients))
  names(primers.after) <- paste0("AfterTrim_", names(primers.after))
  # Combine with quality filter before/after
  qa.out <- cbind(out[1:2,], qaSummary.after.trim, primers.before, primers.after)
  }
  colnames(qa.out)[1:3] <- c("inputReads","filteredReads","trimmedReads")
  qa.out <- qa.out[,colSums(qa.out) != 0]
} else cat("Amplicon must be specified as '16S' or 'ITS.'")
  toc()
  return(qa.out)
}

