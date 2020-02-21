# Function for merging sequence tables from multiple runs

#load ASV tables and merge, set output path.----
#this needs a lot of memory.

mergeRemoveChimeras <- function(seqRuns, output.path, method = "consensus", multithread = TRUE, ...){
  
#Merge multiple runs
seqTabs <- list()
samples <- list()
for (s in 1:length(seqRuns)){
  seqTabs[[s]] <- readRDS(seqRuns[[s]])
  samples[[s]] <- rownames(seqTabs[[s]])
  cat(paste0("\n\nDimensions of ESV table for ", basename(seqRuns[[s]]), ":\nSamples: ", dim(seqTabs[[s]])[1], "\nESVs: ",dim(seqTabs[[s]])[2]))
}

# find duplicated sampleIDs, and append the sequenceRunID to those names.
duped <- unlist(samples)[duplicated(unlist(samples))]
for (s in 1:length(seqRuns)){
 # for (s in 1:3){
    rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped] <-
    paste0(rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped], "_", substr(basename(seqRuns[s]), 1, 5))
    #print(paste0(rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped], "_", substr(basename(seqRuns[s]), 1, 5)))
}

cat("\nMerging ESV tables...")
tic()
st.all <- mergeSequenceTables(tables = seqTabs)
toc()
cat(paste0("\nDimensions of merged ESV table:\nSamples: ", dim(st.all)[1], "\n\nESVs: ",dim(st.all)[2]))

#Remove chimeras
tic()
cat("\nRemoving chimeras...")
seqtab <- st.all
seqtab.nochim <- removeBimeraDenovo(st.all, method=method, multithread=multithread, verbose=TRUE)
toc()
cat(paste0("\nPercent of non-chimeric reads:", sum(seqtab.nochim)/sum(st.all)))
saveRDS(seqtab.nochim, output.path)
cat(paste0("\nMerged ESV table saved to: ",output.path))
return(seqtab.nochim)
}

