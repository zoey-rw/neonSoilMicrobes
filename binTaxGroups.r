# get abundances by taxon rank.
get_tax_level_abun <- function(ps, tax_rank_list = c("phylum","class","order","family","genus"), min_seq_depth = 1000) {
require("phyloseq")
  require("speedyseq")
  require("data.table")
  ps_orig <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
  out.list <- list()
for (r in 1:length(tax_rank_list)) {
  ps <- ps_orig
  tax_rank <- tax_rank_list[[r]]
  cat(paste("\nEvaluating at rank:", tax_rank))
  # get sequence total
  seq_total <- rowSums(otu_table(ps))
  
  if (!tax_rank %in% c("phylum","class","order","family","genus")) {
    tax_table(ps) <-  tax_table(ps)[,tax_rank]
  }
  
  glom <- speedyseq::tax_glom(ps, taxrank=tax_rank) 
glom_melt <- speedyseq::psmelt(glom)
form <- as.formula(paste0("sampleID ~ ", tax_rank))
  
glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
out_abun <- out_abun[sample_names(ps),]
out_abun$other <- NULL
out_abun$other <- seq_total - rowSums(out_abun)

# turn into relative abundances
out_rel <- out_abun/rowSums(out_abun)

# Make a prevalence (frequency) table too
# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(glom), 2, function(x){sum(x > 0)})
N.SVs <- data.frame(table(tax_table(ps)[, tax_rank], exclude = NULL))
colnames(N.SVs)[2] <- "N.SVs"
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(prevalence = prevdf/nsamples(glom),
                    totalAbundance = taxa_sums(glom),
                    tax_table(glom)[,tax_rank])
prevdf <- merge(prevdf, N.SVs, by.x = colnames(prevdf)[3], by.y  = "Var1", all.x=T)

out <- list(out_abun, out_rel, seq_total, prevdf, ps)
names(out) <- c("abundances","rel.abundances","seq.total","prevalence","phyloseq")
out.list[[r]] <- out
}
names(out.list) <- tax_rank_list
return(out.list)
}

