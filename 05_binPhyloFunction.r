rm(list=ls())
library(data.table)
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/common_group_quantification.r')

#set output path.----
output.path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/neonAbundances.rds"
# load NEON SV table as otu file
otu <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/otuTable.rds")
# load NEON taxonomy
tax <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/taxTable.rds")
# load tax-to-function key
tax_fun <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/pecan_gen/reference_data/bacteria_tax_to_function.rds")


# remove leading "k__" in taxonomy.
for (i in 1:ncol(tax)) {
  tax[, i] <- substring(tax[, i], 4)
}



# for everything to be lower case.
tax <- apply(tax,2,tolower)
tax_fun[,1:2] <- apply(tax_fun[,1:2],2,tolower)
colnames(tax_fun) <- tolower(colnames(tax_fun))
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'bacteria',] # removes ~300 counts
otu <- otu[, colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

otu <- otu[rowSums(otu) >= 5000,]

#check that things are properly ordered.----
if(sum(rownames(tax) == colnames(otu)) != ncol(otu)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

tax.rownames.save <- rownames(tax)
#tax[,1:7] <- apply(tax[,1:7],2,tolower)
#colnames(tax)[1:7] <- tolower(colnames(tax)[1:7])
tax$phylum <- as.character(tax$phylum)
table(tax$phylum)
tax[which(tax$phylum==""),]$phylum <- "unknown"
#get each level of taxonomy output.----
of_interest <- colnames(tax)[!colnames(tax) %in% c("kingdom","species")]
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   tax_level = of_interest[i],
                                                  groups = unique(tax[,colnames(tax) == of_interest[i]]),
                                                   ref_filter = F)
}
names(all_taxa_out) <- of_interest

phylum.out <- all_taxa_out$phylum$rel.abundances
jerc <- phylum.out[grep("JERC", rownames(phylum.out)),]

bart <- phylum.out[grep("BART", rownames(phylum.out)),]

bart.from.neon <- read.csv("https://s3.data.neonscience.org/neon-os-data/data-frame/BART_042-O-4.5-4-20160725-GEN-DNA1_16S_20190606T210859.560Z.csv")




stei <- phylum.out[grep("STEI", rownames(phylum.out)),]
ornl.mine <- phylum.out[grep("ORNL", rownames(phylum.out)),]


stei.from.neon <- readRDS("otu_validation/all.rds")

library(dplyr)
stei.phylum <- stei.from.neon %>% 
  group_by(dnaSampleID, phylum) %>% 
  summarise(n = sum(individualCount))

stei.neon.rel <- stei.phylum %>% group_by(dnaSampleID) %>% 
  mutate(per=round(n/sum(n), 4))


stei.neon.rel <- as.data.frame(stei.neon.rel[,!colnames(stei.neon.rel) %in% c("n")])
wide <- tidyr::spread(stei.neon.rel, key = phylum, value = per)
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

#save output.----
saveRDS(all_taxa_out,output.path) 

