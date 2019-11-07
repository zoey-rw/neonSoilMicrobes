# download community composition data from NEON

site = "all"
startdate="2016-01"
enddate="2019-07"
check.size=F
amplicon="16S"
out.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/otu_validation/"


#### RETRIEVE DOWNLOAD URLS AND METADATA FROM NEON ####


sites_wanted <- c("CPER", "STER", "DSNY", "OSBS", "HARV")
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10081.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)
site_date <- cbind(avail$data$siteCodes[1],avail$data$siteCodes[2])


#grab a vector of the urls to data. One per unique site-date combination.
core.urls <- unlist(avail$data$siteCodes$availableDataUrls)

sites_out<- list()
for(i in 1:nrow(site_date)){
  site <- site_date[i,1]
  if (!site %in% sites_wanted) next()
  dates <- unlist(site_date[i,2])
  dates_out <- list()
  for(k in 1:length(dates)){
    date <- dates[k]
    site.date <- paste0(site,'/',date)
    
    #grab DP1.10108.001 data for a particular site-date combination.
    core.JSON  <- httr::GET(core.urls[grep(site.date, core.urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
    
    # subset to only the relevant files
    keep <- core.files$data$files$url[!grepl("expanded|Metadata|variables|EML|readme|validation|basic|ITS",core.files$data$files$url)]
    if (length(keep)==0) next()
    
    core.data <- list()
    for (u in 1:length(keep)){
      core.data[[u]]  <- read.delim(keep[u], sep=",")
     if (nrow(core.data[[u]])==0) next()
     core.data[[u]]$dateID <- date
      core.data[[u]]$siteID <- site
    }
      site_date.output <- plyr::rbind.fill(core.data)
      if(nrow(site_date.output)==0) next()
      dates_out[[k]] <- site_date.output
  }
  site_out <- do.call(rbind, dates_out)
    sites_out[[i]] <- site_out
}

all <- do.call(rbind, sites_out)
saveRDS(all, file = "otu_validation/all.rds")



comp.from.neon <- readRDS("otu_validation/all.rds")

library(dplyr)
neon.phylum <- comp.from.neon %>% 
  group_by(dnaSampleID, phylum) %>% 
  summarise(n = sum(individualCount))

neon.phylum.rel <- neon.phylum %>% group_by(dnaSampleID) %>% 
  mutate(per=round(n/sum(n), 4))


neon.phylum.rel <- as.data.frame(neon.phylum.rel[,!colnames(neon.phylum.rel) %in% c("n")])
wide <- tidyr::spread(neon.phylum.rel, key = phylum, value = per)
wide[is.na(wide)] <- 0
head(as.data.frame(wide))
rowSums(wide[,2:35])

wide$siteID <- substr(wide$dnaSampleID, 1, 4)
wide$plotID <- substr(wide$dnaSampleID, 1, 8)
wide$dnaSampleID <- as.character(wide$dnaSampleID)
wide$dateID <- sapply(strsplit(basename(wide$dnaSampleID), split="-"), '[', 5)

ster <- wide[grep("STER", wide$dnaSampleID),]
osbs <- wide[grep("OSBS", wide$dnaSampleID),]
cper <- wide[grep("CPER", wide$dnaSampleID),]
dsny <- wide[grep("DSNY", wide$dnaSampleID),]
harv <- wide[grep("HARV", wide$dnaSampleID),]
saveRDS(wide, "otu_validation/rel_abundance_5sites.rds")
