# source functions, load libraries
source("helperFunctions.r")
source("dada2.R")
source("trimPrimers.R")
# set parameters for analysis
amplicon <- "16S"
wd <- getwd() # change this path if you want data to download elsewhere

#### CREATE OUTPUT LOCATIONS
# Creates output folders for each step of analysis
create_data_directory(path = wd) # creates a "data" directory with subdirectories
# Output for raw sequences
raw.seqs <- "data/raw_seqs/16S/fastq/"
# set output directory for trimmed sequences (separated by sequencing run)
trimmed.seqs <- "data/trimmed_seqs/16S/"
# set output directory for filtered sequences (separated by sequencing run)
filt.seqs <- "data/filt_seqs/16S/"
# set output directory for sequence tables (separated by sequencing run)
seq.tabs <- "data/seq_tables/16S/"

#### DOWNLOAD RAW SEQUENCES ####
samples <- downloadRawSequenceData(site = "HARV", amplicon = amplicon, 
                        startdate="2017-05", enddate="2017-06", 
                        out.dir = "/data/raw_seqs/ITS", check.size=F)

#### TRIM PRIMERS ####
seqRuns <- sort(list.files(raw.seqs))
track.file <- file.path(trimmed.seqs, paste0("trackPrimers_", amplicon, ".txt")) # tracking file
out <- list()
for(s in 1:length(seqRuns)) {
  raw.seq.dir <- file.path(raw.seqs, seqRuns[s])  # raw file directory for this run
  out.dir <- file.path(trimmed.seqs, seqRuns[s])  # primer-trimmed file directory for this run
  cat(paste("\n\nTrimming primers for sequencing run:", seqRuns[s],"\n\n\n"))
  out[[s]] <- trimPrimers(raw.seq.dir = raw.seq.dir,
                          #samples = samples,
                          out.dir = out.dir)
  names(out)[[s]] <- seqRuns[[s]]
}
all.out <- rbind.named.dfs(out)
write.table(all.out, track.file)

#### QUALITY-FILTER ####
track.file <- file.path(trimmed.seqs, paste0("trackPrimers_", amplicon, ".txt")) # tracking file
out <- list()
for(s in 1:length(seqRuns)) {
  trimmed.seq.dir <- file.path(trimmed.seqs, seqRuns[s])  # trimmed file directory for this run
  out.dir <- file.path(filt.seqs, seqRuns[s])  # quality-filtered file directory for this run
  cat(paste("\n\nTrimming primers for sequencing run:", seqRuns[s],"\n\n\n"))
  out[[s]] <- qualityFilter(raw.seq.dir = raw.seq.dir,
                          #samples = samples,
                          out.dir = out.dir)
}


# RUN DADA2 TO INFER EXACT SEQUENCE VARIANTS (ESVs)
for (s in 1:length(seqRuns)) {
  cat(paste("Running dada2 inference for sequencing run:", seqRuns[[s]]))
  filtpath <- file.path(filt.seqs, seqRuns[[s]])
  out.file <- paste0(out.dir, "/", seqRuns[s],".rds")
  runDada2(filtpath, out.file, seed = 100)
  cat(paste("Exact sequence variant table constructed for ",  seqRuns[[s]],". Saving...\n"))
}
