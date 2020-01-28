# Function to download 16S and ITS raw sequencing data from NEON soil samples. 
#
# Initially written by Lee Stanish, NEON; further adapted by Zoey Werbin
# Will not work for legacy data.
#
#downloadRawSequenceData(out.dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/testing/raw_seqs/ITS", startdate="2016-01", enddate="2016-04")

site = c("CPER","DSNY","HARV","OSBS","STER"); 
startdate="2015-06"; 
enddate="2019-06"; 
out.dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/testing/raw_seqs/ITS"; 
tar_dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/ITS/rawSeqs/recent/tar/";  
amplicon="ITS"; fastq_dir = NULL; check.size=F; 

site = c("CPER","DSNY","HARV","OSBS","STER"); 
startdate="2015-06"; 
enddate="2019-06"; 
out.dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/testing/raw_seqs/16S"; 
tar_dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/16S/rawSeqs/recent/tar/"; 
amplicon="16S"; check.size=F; fastq_dir = NULL;  

downloadRawSequenceData(site = "CPER", amplicon = "16S", startdate="2016-06", enddate="2016-07", check.size=F)

downloadRawSequenceData <- function(site = "HARV", startdate="YYYY-MM", enddate="YYYY-MM", check.size=F, amplicon="ITS", out.dir= NULL, fastq_dir = NULL, tar_dir = NULL, specific.runs=NULL) {
  
  source("helperFunctions.r")
  
  #### CHECK INPUTS ####
  if (!require(neonUtilities)) install.packages('neonUtilities')
  if(amplicon != "16S" & amplicon != "ITS"){
    stop("Value of 'amplicon' must be '16S' or 'ITS'. 16S sequencing targets bacterial and archaeal sequences, while ITS sequencing targets fungal sequences.")
  }
  
  #### CREATE OUTPUT DESTINATIONS ####
  if (is.null(out.dir)) out.dir <- getwd()
  if (!dir.exists(out.dir)) dir.create(out.dir, showWarnings = FALSE)
  if (is.null(fastq_dir)) fastq_dir <- file.path(out.dir, "fastq/")
  if (is.null(tar_dir)) tar_dir <- file.path(out.dir, "tar/")
  if (!dir.exists(fastq_dir))  dir.create(fastq_dir, showWarnings = FALSE)
  if (!dir.exists(tar_dir)) dir.create(tar_dir, showWarnings = FALSE)
  
  #### RETRIEVE DOWNLOAD URLS AND METADATA FROM NEON ####
  
  # Load amplicon sequencing metadata into workspace
  dat <- neonUtilities::loadByProduct(dpID = "DP1.10108.001", site = site, startdate = startdate,
                                      enddate = enddate, package = "expanded", check.size = check.size)
  # Extract soil raw file metadata.
  sequencing.dat <- dat[["mmg_soilRawDataFiles"]]
  # Subset by amplicon type.
  sequencing.dat <- sequencing.dat[grep(amplicon, sequencing.dat$rawDataFileName), ]
  # Add useful columns
  sequencing.dat$which.reads <- substr(sequencing.dat$rawDataFileDescription, 1, 2) # get "R1" or "R2" from column description
  sequencing.dat$new_names <-  paste0(sequencing.dat$dnaSampleID,  "_", amplicon, "_", sequencing.dat$which.reads, ".fastq")
  # Preferentially remove duplicates that don't have an internalLabID
  sequencing.dat <- sequencing.dat[order(sequencing.dat$internalLabID),]
  sequencing.dat <- sequencing.dat[!duplicated(sequencing.dat$new_names),]
  
  # Grab internalLabID from other dataframe if needed
  other.dat <- dat[[paste0("mmg_soilMarkerGeneSequencing_", amplicon)]]
  missing.internalID <- sequencing.dat[is.na(sequencing.dat$internalLabID),]$dnaSampleID
  fix.internalID <- other.dat[other.dat$dnaSampleID %in% missing.internalID, c("dnaSampleID", "internalLabID")]
  colnames(fix.internalID)[2] <- "newinternalID"
  sequencing.dat <- merge(sequencing.dat, fix.internalID, by = "dnaSampleID", all.x=T)
  sequencing.dat$internalLabID <- ifelse(is.na(sequencing.dat$internalLabID), sequencing.dat$newinternalID, sequencing.dat$internalLabID)
  sequencing.dat <- sequencing.dat[!is.na(sequencing.dat$internalLabID ),]
  
  # Create file names to search against our already-downloaded files
  sequencing.dat$fastqFileNames <- paste0(sequencing.dat$internalLabID, "_", amplicon, "_", sequencing.dat$which.reads, ".fastq")
  
  cat(paste0("File metadata successfully downloaded. Approximately ", length(unique(sequencing.dat$fastqFileNames)), " unpaired FASTQ files hosted by NEON.\n"))
  
  
  # Get unique tar download URLs
  u.urls <- sort(unique(sequencing.dat$rawDataFilePath))
  n.urls <- length(u.urls)
  # get tar file names
  fileNms <- gsub('^.*\\/', "", u.urls)
  seqRun <- sapply(strsplit(fileNms, split = "_"), '[', 2)
  prev.downloaded <- length(intersect(fileNms, list.files(tar_dir, pattern = ".tar.gz")))
  if (prev.downloaded > 0) {
    cat(paste0("\n\nDesired FASTQ files are split into ", n.urls, " unique sequencing runs (.tar.gz files), each of which is .5-3 GB. Current `tar_dir`` has ", prev.downloaded, " of ", n.urls, " tar.gz files.")) 
  } else cat("If tar.gz files have previously been downloaded, set the `tar_dir` argument to the previous download location.")
  
  
  #### LOOP THROUGH SEQUENCING RUNS AND DOWNLOAD DESIRED FILES ####
  
  for (i in 1:length(u.urls)) {
    #  for (i in 10:10) {
    if(!is.null(specific.runs)) if ((!seqRun[i] %in% specific.runs) & (!i %in% specific.runs)) next()
    # check if file needs to be downloaded
    file_exists <- grepl(fileNms[i], list.files(tar_dir), fixed = TRUE)
    if (any(file_exists)) {
      cat(paste0("\n\nRun: ", i, "/", n.urls, ": ", fileNms[i], " already exists. Skipping download."))
    } else {
      
      #### DOWNLOAD tar files ####
      utils::download.file(url = as.character(u.urls[i]), destfile = paste0(tar_dir, fileNms[i]))
      cat(paste0("\nRun: ", i, "/", n.urls, ": ", fileNms[i], " has been downloaded to:\n", paste0(tar_dir, fileNms[i])))
    }
    
    # Create output directory for fastQ files from this run
    fastq.out.dir <- paste0(fastq_dir, seqRun[i])
    dir.create(fastq.out.dir, showWarnings = FALSE)
    
    # Fix missing .tar extension 
    if (!grepl("tar", fileNms[i])) {
      base.run <- strsplit(fileNms[i], ".fastq")[[1]][1]
      file.rename(from = paste0(tar_dir, base.run, ".fastq.gz"), 
                  to = paste0(tar_dir, base.run, "_fastq.tar.gz"))
    }
    
    # Get overlap between our tar files and our metadata files
    files.in.tar <-  utils::untar(tarfile = paste0(tar_dir, fileNms[i]), list = T)
    search.list  <-  paste0(sequencing.dat$fastqFileNames, collapse = "|")
    files.to.extract <- files.in.tar[grep(search.list, files.in.tar)]
    
    # Get names of files we've already downloaded; remove from search list
    already_have <- list.files(fastq.out.dir, recursive = T)
    already_have_processed <- already_have[grepl("-GEN-DNA", already_have)] # subset to the *new* files from this run.
    already_have_skip <- sequencing.dat[sequencing.dat$new_names %in% already_have_processed,]$fastqFileNames
    files.to.extract <- files.to.extract[!basename(files.to.extract) %in% already_have_skip]
    if(length(files.to.extract) == 0) { 
      cat(paste0("\nNo new files from run: ", i, "/", n.urls, ": ", fileNms[i]))
      next()
    }
    
    # Extract select files from tar
    utils::untar(
      tarfile = paste0(tar_dir, fileNms[i]),
      files = files.to.extract,
      exdir = fastq.out.dir
    )
    
    cat(paste0("\n", length(files.to.extract)," files from run: ", i, "/", n.urls, ": ", fileNms[i], " have been extracted."))
    
    # Move nested files into a new folder
    cmd <-  paste0("find ", fastq.out.dir, " -mindepth 2 -type f -exec mv -i -f '{}' ", fastq.out.dir, " ';'");  system(cmd)
    # Now remove only the empty directories
    cmd <- paste0("cd ", fastq.out.dir, "; find . -type d -empty -delete"); system(cmd) 
    
    #### FIX FILENAMES (i.e. to SITE_PLOT_HORIZON_CORE_DATE_DNA1_R1.fastq) ####
    
    files_extracted <- basename(files.to.extract)
    samples_extracted <- gsub("_ITS_R[12].fastq|_16S_R[12].fastq", "", files_extracted)
    
    # Subset file metadata to samples we have, and reorder
    file.data <-  sequencing.dat[which(sequencing.dat$fastqFileNames %in% files_extracted), ] # subset file metadata to our DNA samples
    file.data <- file.data[order(match(file.data$fastqFileNames, files_extracted)), ] # match order
    file.data <- file.data[!duplicated(file.data$fastqFileNames),] # some internalLabIDs are repeated, like for BMI_Plate35WellG8
    
    # Create old/new file names and rename
    new.names <- file.path(fastq.out.dir, file.data$new_names) 
    new.names <- new.names[!duplicated(new.names)]
    old.names <- file.path(fastq.out.dir, files_extracted)
    if (nrow(file.data) > 0 & 
        (length(old.names) == length(new.names)) &
        all(file.data$fastqFileNames==files_extracted) &
        all(unlist(lapply(old.names, file.exists)))) { 
      file.rename(from = old.names, to = new.names)
      cat(paste0("\nDownloaded the following files:\n", paste0(new.names, collapse="\n")))
    } else cat(paste0("\nFile renaming not successful for run: ", i, "/", n.urls, ": ", fileNms[i]))
    
    # remove any NA files - shouldn't need this anymore
    # system(paste0("rm ", fastq.out.dir, "/NA"))
  }
  cat("\n\nFile download complete.")
}
