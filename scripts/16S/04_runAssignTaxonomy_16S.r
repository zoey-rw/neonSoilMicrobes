# merge all ESV tables frmo each NEON run
source("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/helperFunctions.r")
library(dada2)
library(foreach)
library(doParallel)
cores=detectCores()
registerDoParallel(cores)

#specify output path here.
tax_output_path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/taxTable.rds"
otu_output_path <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/otuTable.rds"
greengenes.path <- "/projectnb/talbot-lab-data/NEFI_data/gg_13_8_train_set_97.fa"
seq.dir <- "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/seqTab/"

seqRuns <- list.files(seq.dir, full.names = T)
    
    # Merge multiple runs
    # seqTabs <- list()
    # samples <- list()
    # for (s in 1:length(seqRuns)){
    #   print(s)
    #   seqTabs[[s]] <- readRDS(seqRuns[[s]])
    #   samples[[s]] <- rownames(seqTabs[[s]])
    #   print(paste0("Dimensions of ESV table for ", basename(seqRuns[[s]]), ":"))
    #   print(dim(seqTabs[[s]]))
    # }
    # 
    # # find duplicated sampleIDs, and append the sequenceRunID to those names. 
    # duped <- unlist(samples)[duplicated(unlist(samples))]
    # for (s in 1:length(seqRuns)){
    #   rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped] <- 
    #     paste0(rownames(seqTabs[[s]])[rownames(seqTabs[[s]]) %in% duped], "_", substr(basename(seqRuns[s]), 1, 5))
    # }
    # st.all <- mergeSequenceTables(tables = seqTabs)
    # print(paste0("Dimensions of merged ESV table:"))
    # print(dim(st.all))
    # 
    # #Remove chimeras
    # tic()
    # cat("Removing chimeras...")
    # seqtab <- st.all
    # seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
    # sum(seqtab.nochim)/sum(st.all)
    # toc()
    # saveRDS(seqtab.nochim, otu_output_path)
    
    seqtab.nochim <- readRDS(otu_output_path)

    # test.otu <- seqtab.nochim[,1:10000]
    test.to_assign <- colnames(seqtab.nochim)
    #test.to_assign <- colnames(seqtab.nochim)[1:20000]
    sequence.list <- split(test.to_assign, ceiling(seq_along(test.to_assign)/10000))

    tax.list <- list()
    tax.list <- foreach(s=1:length(sequence.list), .errorhandling='pass') %dopar% {
      result <- tryCatch({
        
          #assign taxonomy.
          tic()
          cat('Assigning taxonomy using the greengenes Classifier...\n')
          to_assign <- sequence.list[[s]]
          out <- dada2::assignTaxonomy(to_assign,greengenes.path,multithread = T, verbose = TRUE)
          tax.list[[s]] <- out
          cat('Taxonomy assignment complete. ')
          toc()
          out
            }, error=function(e){paste0("ERROR in assigning taxonomy to chunk #: ", s,"; skipping to next chunk.")})
            return(result) #returned to doParallel
    }
    saveRDS(tax.list, "/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/taxTable_list.rds")
    
    tax <- do.call(rbind, tax.list)
    #save output as your taxonomy file.
    saveRDS(tax, tax_output_path)
    