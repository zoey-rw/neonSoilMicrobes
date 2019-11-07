# Function to download 16S raw sequencing data from NEON soil samples. 
#
# Initially written by Lee Stanish, NEON; further adapted by Zoey Werbin
# Will not work for legacy data.
#
# for testing:
startdate="2016-01"
enddate="2019-07"
check.size=F
site="CPER"
outdir="/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/"
amplicon="16S"

downloadRawSequenceData(site="CPER", startdate="2013-01", enddate="2019-01", check.size=F, outdir="/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/", amplicon = "16S")
downloadRawSequenceData(site="HARV", startdate="2013-01", enddate="2019-01", check.size=F, outdir="/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/", amplicon = "16S")
downloadRawSequenceData(site="BART", startdate="2013-01", enddate="2019-01", check.size=F, outdir="/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/", amplicon = "16S")

downloadRawSequenceData(site="all", startdate="2013-01", enddate="2019-01", check.size=F, outdir="/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/raw_sequences_16S/recent/", amplicon = "16S")

downloadRawSequenceData <- function(site="all", startdate="YYYY-MM", enddate="YYYY-MM", outdir="", check.size=TRUE, 
                                    amplicon = "16S", ...){
  
  if (!require(neonUtilities)) install.packages('neonUtilities')
  if(amplicon != "16S" & amplicon != "ITS"){
    stop("Value of 'amplicon' must be '16S' or 'ITS'. 16S sequencing targets bacterial and archaeal sequences, while ITS sequencing targets fungal sequences. Not currently tested for ITS data.")
  }
  
  #### CREATE OUTPUT DESTINATIONS ####
  fastq_dir <- file.path(outdir,"fastq/")
  tar_dir <- file.path(outdir,"tar/")
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
  cat(paste0("There are ", length(u.urls), " unique files to download.\n")) 
  
  for(i in 1:length(u.urls)) { # loop through each URL
    #for(i in 2:3) { # for testing
    
    # search for the tar file name in the provided out-directory
    file_exists <- grepl(fileNms[i], list.files(tar_dir), fixed = TRUE)
    if(any(file_exists)) { cat(paste0("File: ", fileNms[i], " already exists. Skipping download.\n"))

    } else {
      
      #### DOWNLOAD TAR FILES AND EXTRACT THE SAMPLES WE WANT ####
      utils::download.file(url=as.character(u.urls[i]), destfile = paste0(tar_dir, fileNms[i]))
      cat(paste0("File: ", fileNms[i], " has been downloaded.\n"))
      
    }
        # look into compressed tar files and get the list of FASTQ file names
        files_in_tar <- utils::untar(tarfile = paste0(tar_dir, fileNms[i]), list = T)
        files_we_have <- gsub("\\_R2|\\_R1", "", files_in_tar)
        
        # now subset the FASTQ files to only those which are in our metadata
        files_to_extract <- files_in_tar[which(basename(files_we_have) %in% sequencing.dat$files_we_want)]
        # extract and decompress those FASTQ files into a new "decompressed" subdirectory. 
        # this can take a minute.
        if (length(files_to_extract) > 0) {
          cat("Extracting files...")
          utils::untar(tarfile = paste0(tar_dir, fileNms[i]), 
                files = files_to_extract,
                exdir = paste0(fastq_dir))
          cat(paste0("File: ", fileNms[i], " has been extracted.\n"))
          
          
          #### FIX FILENAMES (i.e. to SITE_PLOT_HORIZON_CORE_DATE_DNA1_R1.fastq) ####
          
          # create sample names with R1 or R2 prefix 
          R1 <- grep("\\_R1", files_to_extract)
          R2 <- grep("\\_R2", files_to_extract)
          sequencing.dat$new_names <- gsub("\\-GEN", "", sequencing.dat$dnaSampleID) # don't need the "GEN" suffix 
          if (length(R1) > 0) {
            sequencing.dat$new_names <- paste0(sequencing.dat$new_names,"_R1.fastq")
          } else if(length(R2) > 0) {
            sequencing.dat$new_names <- paste0(sequencing.dat$new_names,"_R2.fastq")} # this is not a solid way of appending R1/R2
          
          # find the metadata row that matches our downloaded sample files
          new_names <- sequencing.dat[charmatch(substr(basename(files_to_extract), 1, 19), sequencing.dat$processedSeqFileName),]$new_names
          
          # create old/new file names and rename
          old.names <- paste0(outdir, "fastq/", files_to_extract)
          new.names <- paste0(outdir, "fastq/", dirname(files_to_extract), "/", new_names)
          file.rename(from = old.names, to = new.names)
          
          
          #### CLEAN UP TIME! ####  
          
          # moves nested files out of their folders, and into a new folder called "to_keep"
          cmd <- paste0("find ", fastq_dir, " -mindepth 2 -type f -exec mv -i -f '{}' ", fastq_dir, " ';'")
          system(cmd)
          # now let's remove only the empty directories
          cmd <- paste0("rmdir -p ", fastq_dir, unique(dirname(files_to_extract)))
          system(cmd) # ignore the error here! (I don't know how to avoid it)
          
          cat("The following samples have been renamed and are in the directory 'fastq/':")
          print(basename(new.names))
          
        } else {
          cat("No fastq files to extract.")
        } # close conditional extraction of new/wanted fastq files
    } # close loop through URLs
  } ## close function ##
  