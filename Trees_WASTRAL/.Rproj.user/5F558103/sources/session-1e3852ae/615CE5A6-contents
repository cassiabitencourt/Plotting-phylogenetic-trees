# --- Code written by Cassia Bitencourt ---
# This script if for tree topology comparison.

# Install and load required packages
# install.packages("phytools")
# install.packages("phangorn")
# install.packages("ape")
library(phytools)
library(phangorn)
library(ape)
library(treeio)
getwd()
#setwd()

# Load trees with branch lenghts optimised after wASTRAL with no annotation
all_raxml_s0_dt2 <- read.tree("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/dataset2_diamond/all_raxml_AAtoDNA/BrLe_SpeciesTree_PxrrRooted_raxml_allebg_wastral.s0.full.tre.raxml.bestTree")
data <- read.csv("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/names_dt1.csv")
all_raxml_s0_dt2 <- rename_taxa(tree = all_raxml_s0_dt2, data, key = 1, value = 2)

nha_raxml_s0_dt2 <- read.tree("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/dataset2_diamond/nha_raxml_AAtoDNA/BrLe_SpeciesTree_PxrrRooted_dt2_raxml_wastral.s0.nhaebg.tre.raxml.bestTree")
data <- read.csv("/Users/c.bitencourt/Documents/manuscripts/phylogenomics_Apocynaceae/Trees_WASTRAL/Data/names_dt1.csv")
nha_raxml_s0_dt2 <- rename_taxa(tree = nha_raxml_s0_dt2, data, key = 1, value = 2)

# Ladderize the trees before creating the tanglegrams
all_raxml_s0_dt2 <- ladderize(all_raxml_s0_dt2)
nha_raxml_s0_dt2 <- ladderize(nha_raxml_s0_dt2)

# RAxML-adaptive object
# Use the cophylo function to create a tanglegram for IQTrees
cophylo_obj_rx_dt2 <- cophylo(all_raxml_s0_dt2, nha_raxml_s0_dt2)#, assoc_iq_dt2)

# Save the plot to a PDF
pdf("Results/tanglegram_plot_allvsnha_raxml_s0_dt2.pdf", width = 10, height = 35)

# Adjust margins for titles
par(mar = c(4, 4, 7, 4)) # Increase top margin to fit titles

# Plot the tanglegram
plot(cophylo_obj_rx_dt2,
     link.type = "curved",       
     link.lwd = 2,            
     link.col = "darkred",
     alpha = 0.6,
     fsize = 0.15,               # Font size for tip labels
     link.lty = 1,               # Line type for links
     rotate = TRUE,              # Rotate tip labels for better readability 
     cex = 5,                    # Scale of the tree
     offset = 10,                # Adjust offset for tip labels
     link.offset = 0.01,         # Decrease to shorten the link line length
     x.lim = c(5, 3000),         # Extend x-axis limits
     y.lim = c(5, 1500))         # Extend y-axis limits for titles

# Close the PDF device
dev.off()

