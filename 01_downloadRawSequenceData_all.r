# script to download and rename everything from NEON, but keep it separated by sequencing run.

site = "all"
startdate="2016-01"
enddate="2019-07"
check.size=F
amplicon="16S"
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/"


#### CREATE OUTPUT DESTINATIONS ####
fastq_dir <- file.path(out.dir,"fastq/")
#fastq_dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/test/"
tar_dir <- file.path(out.dir,"tar/")
#tar_dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/test/tar/"
if(!dir.exists(fastq_dir))  dir.create(fastq_dir)
if(!dir.exists(tar_dir))    dir.create(tar_dir)

#### RETRIEVE DOWNLOAD URLS AND METADATA FROM NEON ####

# load amplicon sequencing metadata into workspace
dat <- neonUtilities::loadByProduct(dpID="DP1.10108.001", site=site, startdate=startdate, enddate=enddate, package="expanded", check.size=check.size)

# extract 16S or ITS marker gene sequencing metadata.
sequencing.dat <- dat[[paste0("mmg_soilMarkerGeneSequencing_", amplicon)]]

# clean up some FASTQ filenames
sequencing.dat$files_we_want <- gsub(".gz","", sequencing.dat$processedSeqFileName)
sequencing.dat$files_we_want <- gsub("\\_fastq", "\\.fastq", sequencing.dat$files_we_want)

# get unique tar download URLs, and subset by amplicon type.
u.urls <- unique(dat$mmg_soilRawDataFiles$rawDataFilePath)
u.urls <- u.urls[grepl(amplicon, u.urls, fixed = TRUE)]

# get tar file names
fileNms <- gsub('^.*\\/', "", u.urls)
seqRun <- sapply(strsplit(fileNms, split="_"), '[', 2)
cat(paste0("There are ", length(u.urls), " unique files to download.\n")) 

for(i in 1:length(u.urls)) { # loop through each URL
  #for(i in 34:length(u.urls)) { # for testing

  # search for the tar file name in the provided out-directory
  file_exists <- grepl(fileNms[i], list.files(tar_dir), fixed = TRUE)
  if(any(file_exists)) { cat(paste0("File: ", fileNms[i], " already exists. Skipping download.\n"))
    
  } else {
    
    #### DOWNLOAD TAR FILES AND EXTRACT THE SAMPLES WE WANT ####
    utils::download.file(url=as.character(u.urls[i]), destfile = paste0(tar_dir, fileNms[i]))
    cat(paste0("File: ", fileNms[i], " has been downloaded.\n"))
    
  }
    utils::untar(tarfile = paste0(tar_dir, fileNms[i]), 
                 exdir = paste0(fastq_dir, seqRun[i]))
    cat(paste0("File: ", fileNms[i], " has been extracted.\n"))
    
    # get names of files from this sequencing run
    files_extracted <- list.files(paste0(fastq_dir, seqRun[i]), recursive = T)
    files_extracted <- files_extracted[which(nchar(files_extracted) > 55)] # VERY hacky way to subset to the *new* files from this run.
    
    # moves nested files into a new folder 
    cmd <- paste0("find ", paste0(fastq_dir, seqRun[i]), " -mindepth 2 -type f -exec mv -i -f '{}' ", paste0(fastq_dir, seqRun[i]), " ';'")
    system(cmd)
    # now let's remove only the empty directories
    cmd <- paste0("rmdir -p ", paste0(fastq_dir, seqRun[i]), "/", unique(dirname(files_extracted)))
    system(cmd) # ignore the error here! (I don't know how to avoid it)
    
    files_extracted <- basename(files_extracted)

    
    #### FIX FILENAMES (i.e. to SITE_PLOT_HORIZON_CORE_DATE_DNA1_R1.fastq) ####
    
    # create sample names with R1 or R2 prefix 
    R1 <- grep("\\_R1", files_extracted)
    R2 <- grep("\\_R2", files_extracted)
    sequencing.dat$new_names <- sequencing.dat$dnaSampleID 
    if (length(R1) > 0) {
      sequencing.dat$new_names <- paste0(sequencing.dat$new_names,"_R1.fastq")
    } else if(length(R2) > 0) {
      sequencing.dat$new_names <- paste0(sequencing.dat$new_names,"_R2.fastq")} # this is not a solid way of appending R1/R2
    
    # find the metadata row that matches our downloaded sample files
    extractedLabIDs <- sapply(strsplit(files_extracted, split="_16S"), '[', 1)
    
    #unlist(lapply(extractedLabIDs, grep, sequencing.dat$internalLabID, value = TRUE))
    
    new_names <- sequencing.dat[match(extractedLabIDs, sequencing.dat$internalLabID),]$new_names
    # create old/new file names and rename
    old.names <- paste0(fastq_dir, seqRun[i], "/", files_extracted)
    new.names <- paste0(fastq_dir, seqRun[i], "/",  new_names)
    file.rename(from = old.names, to = new.names)
    
    # remove any NA file
    system(paste0("rm ",fastq_dir, seqRun[i], "/NA"))
    
}
    