# get abundances by taxon rank.
get_tax_level_abun <- function(ps, rank.list = c("phylum","class","order","family","genus")) {

  out.list <- list()
for (r in rank.list) {
  rank <- rank.list[r]
  cat(paste("Evaluating at rank:", rank))
  # get sequence total
  seq_total <- rowSums(otu_table(ps))
  glom <- tax_glom(ps, taxrank=rank) 
glom_melt <- psmelt(glom)
form <- as.formula(paste0("sampleID ~ ", rank))
glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance")
out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
out_abun <- out_abun[sample_names(ps),]
out_abun$other <- seq_total - rowSums(out_abun)

# turn into relative abundances
out_rel <- out_abun/rowSums(out_abun)

# Make a prevalence (frequency) table too
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(glom), 2, function(x){sum(x > 0)})
N.SVs <- data.frame(table(tax_table(ps)[, rank], exclude = NULL))
colnames(N.SVs)[2] <- "N.SVs"
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(prevalence = prevdf/nsamples(glom),
                    totalAbundance = taxa_sums(glom),
                    tax_table(glom)[,rank])
prevdf <- merge(prevdf, N.SVs, by.x = colnames(prevdf)[3], by.y  = "Var1", all.x=T)

out <- list(out_abun, out_rel, seq_total, prevdf)
names(out) <- c("abundances","rel.abundances","seq.total","prevalence")
out.list[[rank]] <- out
}
return(out.list)
}


test <- get_tax_level_abun(ps)
saveRDS(test, "ps_groupAbundances_ITS.rds")


# compare to colin's function
# source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/common_group_quantification.r')
# 
# #get each level of taxonomy output.----
# of_interest <- colnames(tax)[!colnames(tax) %in% c("kingdom", "species")]
# of_interest <- c("phylum","class")
# # Extract abundance matrix from the phyloseq object
# otu = as.data.frame(as(otu_table(ps), "matrix"))
# tax = as.data.frame(as(tax_table(ps), "matrix"))
# all_taxa_out <- list()
# for(i in 1:length(of_interest)){
#   all_taxa_out[[i]] <- common_group_quantification(otu,
#                                                    tax,
#                                                    tax_level = of_interest[i],
#                                                    groups = unique(tax[,colnames(tax) == of_interest[i]]),
#                                                    ref_filter = F,
#                                                    samp_freq = 0.0)
# }
# names(all_taxa_out) <- of_interest
# 
# colins <- as.data.frame(all_taxa_out$phylum$rel.abundances[order(rownames(all_taxa_out$phylum$rel.abundances)),])
# zoeys <- test$phylum$rel.abundances[order(rownames(test$phylum$rel.abundances)),]
# # Matches exactly
# plot(zoeys$basidiomycota, colins$basidiomycota)
# plot(zoeys$mortierellomycota, colins$mortierellomycota)
# plot(zoeys$ascomycota, colins$ascomycota)
