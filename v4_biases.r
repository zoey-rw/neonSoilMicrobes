rm(list=ls())
library(runjags)
library(foreach)
library(doParallel)
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/tic_toc.r')
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/hierarch_ddirch_means.r')
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/hierarch_ddirch_means.r')
#register parallel environment.----
n.cores <- detectCores()
n.cores <- 16
registerDoParallel(n.cores)

group_abun_16S <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_groupAbundances_16S.rds")


#loop over levels.----
output <- list()
output <-
  foreach(i = 1:length(d)) %dopar% {
    #Get y multivariate matrix.
    abundances <- d[[i]]$abundances
    y <- as.matrix((abundances + 1) / rowSums(abundances + 1)) # need to change this
    
    #get core_plot and plot_site indexing.
    core_plot <- substr(rownames(y), 1, 8)
    core_site <- substr(rownames(y), 1, 4)
    plot_site <- unique(core_plot)
    plot_site <- substr(plot_site, 1, 4)
    
    #fit the hierarchical means.
    fit <- hierarch_ddirch_means(y=y, core_plot = core_plot, plot_site = plot_site, jags.method = 'parallel')
    
    #add row and column names - plot level (you should put this in the function).
    for(j in 1:length(fit$plot.fit)){
      rownames(fit$plot.fit[[j]]) <- unique(core_plot)
      rownames(fit$site.fit[[j]]) <- unique(plot_site)
      colnames(fit$plot.fit[[j]]) <- colnames(y)
      colnames(fit$site.fit[[j]]) <- colnames(y)
    }
    fit$core.fit <- y #add in the core-level data!
    return(fit)
  }
names(output) <- names(d)

#save matrix lists.----
saveRDS(output, output.path)
cat('Script complete. ');
toc()



phys <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10086.soil_phys.rds")
phys <- phys[phys$siteID %in% c("HARV","OSBS","CPER","DSNY","STER"),]
  
ps <- group_abun_16S$ps_16S
samp.df <- sample_data(ps)
phys$soilSampleID <- phys$sampleID
phys$sampleID <- NULL
samp.df <- as.data.frame(samp.df)
df <- as(sample_data(ps), "data.frame")
rownames(df) <- sample_names(ps)
df$soilSampleID <- df$sample
df$sampleID <- NULL
df$dateID <- NULL
df$row.names  <- rownames(df)
phys <- phys[phys$soilSampleID %in% df$soilSampleID,]
samp.df <- merge(df, phys, all.x = T, all.y = F, by = "soilSampleID")
samp.df <- samp.df[!duplicated(samp.df$row.names),]
rownames(samp.df) <- samp.df$row.names
samp.df$sampleID <- samp.df$row.names
samp.df$legacy <- ifelse(samp.df$asDate.x < "2015-01-01", T, F)
sample_data(ps) <- samp.df



abun_phylum <- get_tax_level_abun(ps, tax_rank_list = c("phylum"))
df.plot <- as(sample_data(ps), "data.frame")
df.plot$acido <- abun_phylum$phylum$rel.abundances$acidobacteria
df.plot$color <- ifelse(df.plot$legacy, "red", "blue")
df.plot$source <- ifelse(df.plot$legacy, "NEON_Legacy_V4", "NEON_Recent_V3V4")

df.plot <- df.plot[,c("acido","soilInWaterpH","color","source")]
colnames(df.plot) <- c("acidobacteria","pH","color","source")
#plot(acidobacteria ~ pH, data = df.plot, col = df.plot$color)




#load 16S meta data file with lat/long
d <- readRDS(ramirez_raw_mapping_and_abundance.path)
metadata.cols <- colnames(d)[!grepl("^[dpcofgs]__[bacteria]", colnames(d))]
metadata.cols <- metadata.cols[!grepl("^[dpcofgs]__[vu]", metadata.cols)]
d <- d[,metadata.cols]
d <- d[!d$dataset %in% c("X9","X45","X8"),]
ram.abun <- readRDS(delgado_ramirez_abun.path)
ram.plot <- ram.abun$phylum
ram.plot$sampleID <- rownames(ram.plot)
ram.plot <- merge(ram.plot, d)
ram.plot <- ram.plot[,c("acidobacteria","ph","seq_reg_rev")]
ram.plot$color <- ifelse(ram.plot$seq_reg_rev=="v4", "red", "black")
plot(acidobacteria ~ ph, data = ram.plot, col = ram.plot$color, pch=16)

ram.plot.plot <- ram.plot[,c("acidobacteria","ph","color")]
colnames(ram.plot.plot) <- c("acidobacteria","pH","color")
ram.plot.plot$color <- "green"
ram.plot.plot$source <- "Ramirez"
  
to.plot <- rbind(df.plot, ram.plot.plot)
ggplot(to.plot) + geom_point(aes(x = pH, y = acidobacteria, color = source)) + geom_smooth(aes(x = pH, y = acidobacteria, color = source), span = .5) + ggtitle("Acidobacteria ~ pH, by source") 
#plot(acidobacteria ~ pH, data = to.plot, col = to.plot$color, pch=16, alpha = .3)

colnames(ram.plot) <- c("acidobacteria","pH","source","color")
to.plot <- rbind(df.plot, ram.plot)
to.plot <- to.plot[!to.plot$source %in% c("v1_v2","v5_v7","v3_v5","v1_v3"),]
to.plot$source <- as.factor(to.plot$source)
to.plot <- to.plot[!to.plot$source %in% c("NEON_Recent_V3V4"),]
cbPalette <- c("#000000", "#E69F00", "#56B4E9")

ggplot(to.plot) + geom_point(aes(x = pH, y = acidobacteria, color = source)) + geom_smooth(aes(x = pH, y = acidobacteria, color = source), span = .5) + ggtitle("Acidobacteria ~ pH, by 16S subregion")  + scale_colour_manual(values=cbPalette) + theme_minimal()

# with 4 colors
cbPalette <- c("#000000", "#009E73","#E69F00","#56B4E9")
ggplot(to.plot) + geom_point(aes(x = pH, y = acidobacteria, color = source)) + geom_smooth(aes(x = pH, y = acidobacteria, color = source), span = .5) + ggtitle("Acidobacteria ~ pH, by 16S subregion")  + scale_colour_manual(values=cbPalette) + theme_minimal()



# download marker gene metadata
dat <- neonUtilities::loadByProduct(dpID = "DP1.10108.001", site = c("HARV","OSBS","CPER","DSNY","STER"), startdate = "2018-01",
                                    enddate = "2018-12", package = "expanded", check.size = T)
pcr_old <- dat$mmg_soilPcrAmplification_16S
v4_samples <- pcr_old[pcr_old$targetSubfragment=="v4",]
v3v4_samples <- pcr_old[pcr_old$targetSubfragment=="v3v4",]


# load in all new metadata  - which is v3v4
metadata <- readRDS("data/data_old/allNEONmetadata.rds")
pcr <- metadata$mmg_soilPcrAmplification_16S
new_v3_samples <- pcr[pcr$targetSubfragment=="v3v4",]






map_emp <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/old/prior_abundance_mapping/EMP/emp_map_clean.rds") #environmental metadata
data.list <- readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/old/prior_abundance_mapping/EMP/emp_phylo.level.list_50.path") #relative abundance 

phy_emp <- as.data.frame(t(data.list$phylum))
phy_emp$sampleID <- rownames(phy_emp)
emp_plot <- merge(map_emp, phy_emp)
emp_plot.plot <- emp_plot[,c("Acidobacteria","ph")]
colnames(emp_plot.plot) <- c("acidobacteria","pH")
emp_plot.plot$source <- "EMP_V4"
emp_plot.plot$color <- "black"

to.plot <- rbindlist(list(df.plot, ram.plot.plot, emp_plot.plot),use.names=TRUE)

ggplot(to.plot) + geom_point(aes(x = pH, y = acidobacteria, color = source)) + geom_smooth(aes(x = pH, y = acidobacteria, color = source), span = .5) + ggtitle("Acidobacteria ~ pH, by 16S subregion") 
