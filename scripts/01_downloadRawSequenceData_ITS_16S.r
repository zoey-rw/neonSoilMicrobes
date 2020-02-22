# script to download and rename everything from NEON, but keep it separated by sequencing run.

##### NOTE - after downloading, all the files should be compressed for downstream code to work
# i.e. in the FASTQ directory, type "parallel gzip ::: *" (if GNU parallel is loaded) or "gzip *" otherwise

source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/downloadRawSequenceData.r")

downloadRawSequenceData(site = c("CPER","DSNY","HARV","OSBS","STER"), amplicon = "ITS", startdate="2015-06", enddate="2019-06", out.dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/testing/raw_seqs/ITS", check.size=F, tar_dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data_old/ITS/rawSeqs/recent/tar/")

downloadRawSequenceData(site = c("CPER","DSNY","HARV","OSBS","STER"), amplicon = "16S", startdate="2015-06", enddate="2019-06", out.dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S", check.size=F, tar_dir = "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data_old/16S/rawSeqs/recent/tar/")

    