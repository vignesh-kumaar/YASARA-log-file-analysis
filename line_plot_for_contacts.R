rm(list=ls())
setwd("/home/vignesh/Data/old_ubuntu/YASARA-log-file-analysis")

# filepath <- readline(prompt = "Enter absolute path to the contacts csv file / relative path from working directory ")
# data <- read.csv(filepath, sep="\t", encoding="UTF-8")

library(ggplot2)
library(dplyr)
library(reshape)

line_plot_table <- read.csv("input_files/line_plot_tables/contacts_of_pdb_2.csv", sep="\t")
head(data)


line_plot_table <- line_plot_table %>%
  mutate(X0 = make.unique(as.character(X0)))


line_plot_table$X0 <- factor(line_plot_table$X0, levels = line_plot_table$X0)


ggplot(line_plot_table, aes(x = X0, y = X1, group = 1)) +
  geom_point( alpha = 0.5) +
  labs(x = "Receptor + Binding Residues", y = "Contacts", title = "Line Plot of Binding pairs vs Contacts") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=3))




rm(list=ls())
# Convert the data to a matrix if needed
heat_map_data <- read.csv('input_files/contacts_heat_matrix.csv', sep='\t', header=TRUE)
nrows <- nrow(heat_map_data)
ncols <- ncol(heat_map_data)
heat_map_matrix <- matrix(nrow = nrows, ncol = ncols)
for (i in 1:nrows) {
  for (j in 1:ncols) {
    heat_map_matrix[i, j] <- heat_map_data[i, j]
  }
}

library(lattice)
pal <- colorRampPalette(c("black", "green"), space = "rgb")
levelplot(heat_map_matrix, main="PDB file vs Residue position", xlab="Residue position", ylab="PDB file number", col.regions=pal(8), cuts=8, at=seq(0,70,10), aspect = "fill")



# Comparing interactions_table of two datasets:
rm(list=ls())

# initializing variables:
ecoli_data <- read.csv("ecoli_bamA_mosselii_llpA1/analysis/interactions_table.csv", sep="\t")
pseudomonas_fluorescens_data <- read.csv("input_files/interactions_table.csv", sep="\t")

ecoli_contacts <- as.numeric(unlist(ecoli_data['Contacts.across.PDBs']))
pseudomonas_fluorescens_contacts <- as.numeric(unlist(pseudomonas_fluorescens_data['Contacts.across.PDBs']))

ecoli_hbonds <- as.numeric(unlist(ecoli_data['h_bonds']))
pseudomonas_fluorescens_hbonds <- as.numeric(unlist(pseudomonas_fluorescens_data['h_bonds']))

ecoli_hydrophobic <- as.numeric(unlist(ecoli_data['Hydrophobic']))
pseudomonas_fluorescens_hydrophobic <- as.numeric(unlist(pseudomonas_fluorescens_data['Hydrophobic']))

ecoli_ionic <- as.numeric(unlist(ecoli_data['Ionic']))
pseudomonas_fluorescens_ionic <- as.numeric(unlist(pseudomonas_fluorescens_data['Ionic']))

ecoli_cationpi <- as.numeric(unlist(ecoli_data['CationPi']))
pseudomonas_fluorescens_cationpi <- as.numeric(unlist(pseudomonas_fluorescens_data['CationPi']))

ecoli_pipi <- as.numeric(unlist(ecoli_data['PiPi']))
pseudomonas_fluorescens_pipi <- as.numeric(unlist(pseudomonas_fluorescens_data['PiPi']))


# chisq.test(rbind(pseudomonas_fluorescens_contacts['Contacts.across.PDBs'], ecoli_contacts['Contacts.across.PDBs']))

shapiro.test(pseudomonas_fluorescens_contacts)
shapiro.test(ecoli_contacts)
shapiro.test(pseudomonas_fluorescens_hbonds)
shapiro.test(ecoli_hbonds)
shapiro.test(pseudomonas_fluorescens_hydrophobic)
shapiro.test(ecoli_hydrophobic)
shapiro.test(pseudomonas_fluorescens_ionic)
shapiro.test(ecoli_ionic)
shapiro.test(pseudomonas_fluorescens_cationpi)
shapiro.test(ecoli_cationpi)
shapiro.test(pseudomonas_fluorescens_pipi)
shapiro.test(ecoli_pipi)

# Since data is not normally distributed:
wilcox.test(pseudomonas_fluorescens_contacts, ecoli_contacts)
# p-value reported is 0.3255 > 0.05
# null hypothesis is true, the two data sets are similar

wilcox.test(pseudomonas_fluorescens_hbonds, ecoli_hbonds)

wilcox.test(pseudomonas_fluorescens_hydrophobic, ecoli_hydrophobic)

wilcox.test(pseudomonas_fluorescens_ionic, ecoli_ionic)

wilcox.test(pseudomonas_fluorescens_cationpi, ecoli_cationpi)

wilcox.test(pseudomonas_fluorescens_pipi, ecoli_pipi)


mean(pseudomonas_fluorescens_contacts)
mean(ecoli_contacts)

mean(pseudomonas_fluorescens_hbonds)
mean(ecoli_hbonds)

mean(pseudomonas_fluorescens_hydrophobic)
mean(ecoli_hydrophobic)

mean(pseudomonas_fluorescens_ionic)
mean(ecoli_ionic)

mean(pseudomonas_fluorescens_cationpi)
mean(ecoli_cationpi)

mean(pseudomonas_fluorescens_pipi)
mean(ecoli_pipi)


#
# When lists don't include zeros:
pseudomonas_fluorescens_contacts <- pseudomonas_fluorescens_contacts[pseudomonas_fluorescens_contacts != 0]
ecoli_contacts <- ecoli_contacts[ecoli_contacts != 0]

wilcox.test(pseudomonas_fluorescens_contacts, ecoli_contacts)
# p-value reported is 0.9588 > 0.05

# Permutation test
combined <- c(pseudomonas_fluorescens_contacts, ecoli_contacts)
group <- c(rep("pseudomonas_fluorescens", length(pseudomonas_fluorescens_contacts)), rep("ecoli", length(ecoli_contacts)))

observed_diff <- sum(pseudomonas_fluorescens_contacts) - sum(ecoli_contacts)

perm_diffs <- replicate(10000, {
  perm_group <- sample(group)
  perm_pseduomonas_fluorescens <- combined[perm_group == "pseudomonas_fluorescens"]
  perm_ecoli <- combined[perm_group == "ecoli"]
  sum(perm_pseduomonas_fluorescens) - sum(perm_ecoli)
})

p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
print(p_value)



