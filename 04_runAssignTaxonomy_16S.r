# merge all ESV tables frmo each NEON run

library(dada2)
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/helperFunctions.r")

#specify output path here.
tax_output_path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/taxTable.rds"
otu_output_path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/otuTable.rds"
greengenes.path <- "/projectnb/talbot-lab-data/NEFI_data/gg_13_8_train_set_97.fa"
seq.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_data_construction/seqTab/"

seqRuns <- list.files(seq.dir, full.names = T)

# Merge multiple runs (if necessary)
seqTabs <- list()
samples <- list()
for (s in 1:length(seqRuns)){
  seqTabs[[s]] <- readRDS(seqRuns[[s]])
  samples[[s]] <- rownames(seqTabs[[s]])
  print(paste0("Dimensions of ESV table for ", basename(seqRuns[[s]]), ":"))
  print(dim(seqTabs[[s]]))
}

# find duplicated sampleIDs, and append the sequenceRunID to those names. 
duped <- unlist(samples)[duplicated(unlist(samples))]
for (s in 1:length(seqRuns)){
  rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped] <- 
    paste0(rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped], "_", substr(basename(seqRuns[s]), 1, 5))
}

st.all <- mergeSequenceTables(tables = seqTabs)
print(paste0("Dimensions of merged ESV table:"))
print(dim(st.all))

# Remove chimeras
tic()
cat("Removing chimeras...")
seqtab <- st.all
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(st.all)
toc()

saveRDS(seqtab.nochim, otu_output_path)

#assign taxonomy.
tic()
cat('Assigning taxonomy using the greengenes Classifier...\n')
to_assign <- colnames(seqtab.nochim) #grab sequences to assign.
out <- dada2::assignTaxonomy(to_assign,greengenes.path,multithread = T)
cat('Taxonomy assignment complete. ')
toc()

#save output as your taxonomy file.
saveRDS(out, tax_output_path)


dat$mmg_soilRawDataFiles[which(dat$mmg_soilRawDataFiles$dnaSampleID == "NIWO_004-O-30.5-19.5-20160627-GEN-DNA1"),]
sequencing.dat[which(sequencing.dat$dnaSampleID == "NIWO_004-O-30.5-19.5-20160627-GEN-DNA1"),]