# create taxonomy to function reference dataset.
# currently uses a mix of tidyr, data.table, and base R. 
# 
# This function makes scientifically relevant decisions - edit as needed
# for N-cyclers, taxa that can carry out any entire step of denitrification/nitrification 
# are placed in that group, but taxa with only *some* of the requisite genes are not.
# 
# Also, everything is assigned at the genus level
#
# ref.path is the location of the CSV downloaded from GoogleSheets (lit review by Zoey Werbin, Fall 2018)
# N.path is the location of the CSV shared by Albright
# C.path is the location of the CSV downloaded from the Berlemont supplement

createTaxFunction <- function(ref.path = "data/reference_data/bacteria_func_groups.csv", 
                                     N.path = "data/reference_data/Npathways_Albright2018.csv", 
                                     C.path = "data/reference_data/cellulolytic_Berlemont.csv",
                                     out.path = "data/tax_function_ref.csv"){
require(data.table)
require(tidyr)
require(stringr)
  
# Read in literature-review classifications
fg <- read.csv(ref.path, stringsAsFactors = F)
# Read in Albright et al. N-cycle pathway presence/absence
N_cyclers_raw <-  read.csv(N.path, stringsAsFactors = F)
# Read in Berlemont and Martiny cellulolytic pathway presence/absence
cellulolytic_raw <-  read.csv(C.path, stringsAsFactors = F)

#### 1. Format N-cycle dataset from Albright 2008 ####

# Remove all columns except for taxonomy, environment, and pathways
colnames(N_cyclers_raw)[[1]] <- "Samplename"
N_cyclers <- N_cyclers_raw[,!colnames(N_cyclers_raw) %in% 
                             c("Samplename", "Genome", "StudyName", "Ecosystem", 
                               "Ecosystem.Category", "Ecosystem.Subtype", "Ecosystem.Type",
                               "Environment", "Genome_size_assembled",  "Gene_Count_assembled")]
# Rename some pathways
setnames(N_cyclers,
         c("Nitrogen.Fixation", "Assimilatory.Nitrite.to.ammonia",
           "Dissimilatory.Nitrite.to.Ammonia", "Assimilatory.Nitrate.to.Nitrite",
           "Dissimilatory.Nitrate.to.Nitrite"),
         c("N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction",
           "Assim_nitrate_reduction", "Dissim_nitrate_reduction"))

# Treat "incomplete" pathways as if they are absent
N_cyclers[N_cyclers == "complete"] <- 1
N_cyclers[N_cyclers == "incomplete" | N_cyclers == "None"] <- 0

# Grouping partial nitrification pathway with nitrification,
# and partial denitrification with denitrification.
N_cyclers[N_cyclers$Partial_Nitrification == 1,]$Nitrification <- 1
N_cyclers[N_cyclers$Partial_NO == 1 | N_cyclers$Partial_N2O == 1 |
            N_cyclers$Partial_N2 == 1,]$Denitrification <- 1

# Now we can remove the "partial" columns.
N_cyclers[, c("Partial_Nitrification",
              "Partial_NO",
              "Partial_N2O",
              "Partial_N2")] <- NULL

# Keep only genus-level
N_cyclers$Taxonomic.level <- "Genus"
N_cyclers$Taxon <- N_cyclers$Genus
N_cyclers <- N_cyclers[,!colnames(N_cyclers) %in% c("Domain","Phylum","Class","Order","Family","Genus","Species")]


#### 2. Format dataset of cellulolytic taxa from Berlemont et al. ####

# Subset to genes that are primarily associated with exo- and endo-cellulases (Table 1 from Berlemont paper)
# On the SCC the "Strain" column reads in as "X...Strain" (?)
cellulolytic <- cellulolytic_raw[,colnames(cellulolytic_raw) %in% 
                                   c("Strain","X...Strain", "GH5", "GH6", "GH8", "GH9", "GH12", "GH44", "GH45", "GH48")]
colnames(cellulolytic)[[1]] <- "Strain"
# Create genus column
cellulolytic$genus <- stringr::word(cellulolytic$Strain, 1)
cellulolytic[cellulolytic$genus == "Candidatus",]$genus <- # if first word is 'Candidatus,' grab two words
  stringr::word(cellulolytic$Strain[cellulolytic$genus == "Candidatus"], 1, 2) 

# assign taxa to cellulolytic group
cellulolytic$Cellulolytic <- 0
cellulolytic[apply(cellulolytic, 1, function(x)
  any(x == 1)),]$Cellulolytic <- 1 # if any pathway is present, taxon is cellulolytic
cellulolytic$Taxonomic.level <- "Genus"
cellulolytic$Taxon <- cellulolytic$genus
cellulolytic <- cellulolytic[,c("Taxonomic.level","Taxon","Cellulolytic")]
cellulolytic <- cellulolytic[cellulolytic$Cellulolytic==1,] # drop non-cellulolytic taxa


#### 3. Format groups from literature review ####
fg_lit <- fg[,!colnames(fg) %in% c("Classification.system", "Source", "Notes")]

# Specify which groups to check the spreadsheet for
groups <- c("Nitrification", "Denitrification", "N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction", "Assim_nitrate_reduction", "Dissim_nitrate_reduction","Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph","Copiotroph","Oligotroph")

# Loop through spreadsheet to check for each group
for (g in 1:length(groups)) {
  group <- groups[[g]]
  fg_lit <- fg_lit %>% 
    dplyr::mutate(!!group := dplyr::case_when(Classification == !!group ~ 1,
                                Classification != !!group ~ 0))
}
fg_lit$Classification <- NULL

#### 4. Combine ####
fg_out <- plyr::rbind.fill(cellulolytic, N_cyclers, fg_lit)
fg_out[is.na(fg_out)] <- 0

fg_out <- fg_out[!duplicated(fg_out),]

fg_out[,1:2] <- apply(fg_out[,1:2],2,tolower)
colnames(fg_out) <- tolower(colnames(fg_out))

return(fg_out)
}

test <- createTaxFunction()



