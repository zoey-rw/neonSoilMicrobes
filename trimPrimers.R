# trim primers using cutadapt and the dada2 tutorial instructions.
rm(list=ls())

library(dada2)
library(ShortRead)
library(Biostrings)
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

trimPrimers <- function(raw.seqs, out.dir, flip_rev_primer = FALSE, FWD = "CCTACGGGNBGCASCAG", REV = "GACTACNVGGGTATCTAATCC", cutadapt.path = "/share/pkg.7/cutadapt/1.18/install/bin/cutadapt", test = FALSE, remove.filtN = T){

  seqRun <- list.files(raw.seqs)
  
  source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/helperFunctions.r")
  
# create output directory
if(!dir.exists(out.dir)) dir.create(out.dir)

# read in raw sequences
fnFs <- sort(list.files(raw.seqs, pattern = "_R1.fastq", full.names = TRUE, recursive = T))
fnRs <- sort(list.files(raw.seqs, pattern = "_R2.fastq", full.names = TRUE, recursive = T))

if (test==TRUE){
fnFs <- head(fnFs, 10)
fnRs <- head(fnRs, 10)
}

# remove samples without forward and reverse reads 
fnFs.names <- sapply(strsplit(basename(fnFs), split="_R"), '[', 1)
fnRs.names <- sapply(strsplit(basename(fnRs), split="_R"), '[', 1)
overlap <- intersect(fnFs.names, fnRs.names)
fnFs <- fnFs[fnFs.names %in% overlap]
fnRs <- fnRs[fnRs.names %in% overlap]

# fix mis-specified primer
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

if (flip_rev_primer == T){
REV <- REV.orients[["RevComp"]]
REV.orients <- allOrients(REV)
}

# remove any N reads (dada2 can't handle)
fnFs.filtN <- file.path(out.dir, "filtN", paste0(basename(fnFs), ".gz")) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(out.dir, "filtN", paste0(basename(fnRs), ".gz"))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, truncQ = 2, multithread = TRUE, matchIDs = TRUE)

if (!is.null(cutadapt.path)) system(paste0("export PATH=",cutadapt.path,":$PATH"))
if (system("which cutadapt") == 1) stop("Cannot find path to cutadapt. Set using the 'cutadapt.path' argument, or run from another location.")

# create output file names
fnFs.cut <- file.path(out.dir, basename(fnFs))
fnRs.cut <- file.path(out.dir, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in 1:length(fnFs)) {
  system(paste("cutadapt", R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
               "-m 50",
               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

if(remove.filtN == TRUE){
cat("Removing temporary directory of N-filtered reads...")
  system(paste0("rm -r ", file.path(out.dir, "filtN")))
}

msg1 <- "Primer count before trimming (for first sample):"
primers.before <- checkPrimers(fnFs[[1]], fnRs[[1]], FWD.orients, REV.orients)
qaSummary.before <- qa(fnFs, type="fastq")

msg2 <- "Primer count after trimming (for first sample):"
primers.after <- checkPrimers(fnFs.cut[[1]], fnRs.cut[[1]], FWD.orients, REV.orients)
qaSummary.after <- qa(fnFs.cut, type="fastq")

qa.out <- cbind(head(qaSummary.before[["readCounts"]][,1, drop=F]), head(qaSummary.after[["readCounts"]][,1, drop=F]))
colnames(qa.out) <- c("withPrimersReadCount","primersTrimmedReadCount")
print(qa.out)

return(list(msg1, primers.before, msg2, primers.after, qa.out))
}

