library(data.table)
source("helperFunctions.r")
source("binTaxGroups.r")
source("addBacterialFunction.r")
source("createTaxFunction.r")
library(ggplot2)
abun_ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/ps_groupAbundances_ITS.rds")

plot.list <- list()
for (group in names(abun_ITS)[7:12]){
  samp.df <- sample_data(abun_ITS$ps_ITS)
  y <- abun_ITS[[group]]$rel.abundances[,group]
  samp.df <- samp.df[samp.df$sampleID %in% rownames(abun_ITS[[group]]$rel.abundances),]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  df <- df[df$asDate < "2017-01-01",]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .7) + ggtitle(group) + ylim(0, 1)
}
gridExtra::grid.arrange(plot.list[[1]] + ylim(0, .3), plot.list[[2]], plot.list[[3]], plot.list[[4]])
gridExtra::grid.arrange(plot.list[[5]], plot.list[[6]])

# plot top 8 fungal phyla
cosmo <- abun_ITS$phylum$prevalence[order(-abun_ITS$phylum$prevalence$prevalence),][1:6,c("phylum")]
plot.list <- list()
for (group in cosmo){
  samp.df <- sample_data(abun_ITS$ps_ITS)
  y <- abun_ITS$phylum$rel.abundances[,group]
  samp.df <- samp.df[samp.df$sampleID %in% rownames(abun_ITS$phylum$rel.abundances),]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  df <- df[df$asDate < "2017-01-01",]
  #df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) #+ 
    #theme_classic()
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]],plot.list[[5]], plot.list[[6]])


abun_16S <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/16S/ps_groupAbundances_16S.rds")


plot.list <- list()
for (group in names(abun_16S)[7:19]){
  samp.df <- sample_data(abun_16S$ps_16S)
  y <- abun_16S[[group]]$rel.abundances[,group]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  df <- df[df$asDate < "2017-01-01",]
  #df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) #+ theme_classic()
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])
gridExtra::grid.arrange(plot.list[[5]], plot.list[[6]], plot.list[[7]], plot.list[[8]])
gridExtra::grid.arrange(plot.list[[9]], plot.list[[10]], plot.list[[11]], plot.list[[12]])
plot.list[[13]]




# plot top 10 bacterial phyla
cosmo <- abun_16S$phylum$prevalence[order(-abun_16S$phylum$prevalence$prevalence),][1:10,c("phylum")]
plot.list <- list()
for (group in cosmo){
  samp.df <- sample_data(abun_16S$ps_16S)
  y <- abun_16S$phylum$rel.abundances[,group]
  df <- cbind.data.frame(asDate = samp.df$asDate,
                         siteID = samp.df$siteID,
                         horizon = samp.df$horizon,
                         y = y)
  df <- df[df$asDate < "2017-01-01",]
  #df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
  plot.list[[group]] <- ggplot(df, aes(x = asDate, y = y)) + 
    geom_point(aes(color = siteID)) +
    geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) #+ 
    #theme_classic()
}
gridExtra::grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])
gridExtra::grid.arrange(plot.list[[5]], plot.list[[6]], plot.list[[7]], plot.list[[8]], plot.list[[9]], plot.list[[10]])









# visualize plot means
group <- "oligotroph"
samp.df <- sample_data(abun_16S$ps_16S)
y <- abun_16S[[group]]$rel.abundances[,group]
df <- cbind.data.frame(asDate = samp.df$asDate,
                       siteID = samp.df$siteID,
                       plotID = samp.df$plotID,
                       horizon = samp.df$horizon,
                       y = y)
df <- df[df$asDate < "2017-01-01",]
df <- df %>% group_by(siteID, asDate, plotID) %>% summarize(plot_means = mean(y))
#df <- df[df$asDate < as.Date("20170101", "%Y%m%d"),]
plot.list[[group]] <- ggplot(df, aes(x = asDate, y = plot_means)) + 
  geom_point(aes(color = siteID)) +
  geom_smooth(aes(color = siteID), span = .5) + ggtitle(group) + theme_classic()
