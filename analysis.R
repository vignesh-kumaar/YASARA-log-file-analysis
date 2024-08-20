rm(list=ls())
setwd("/home/dkscara/dkscara/summer-research")

# filepath <- readline(prompt = "Enter absolute path to the contacts csv file / relative path from working directory ")
# data <- read.csv(filepath, sep="\t", encoding="UTF-8")


contacts_table <- read.csv("contacts.csv", sep="\t")
head(data)

library(ggplot2)
library(dplyr)
library(reshape)

residues_and_contacts <- select(contacts_table, Receptor.residue.and.number, Strength.of.contacts)
plot_parameter <- melt(residues_and_contacts, id.vars = "Receptor.residue.and.number")
ggplot(data=contacts_table, aes(y=Strength.of.contacts)) + geom_line()

# ggplot(plot_parameter, aes(x = Receptor.residue.and.number, y = Strength.of.contacts, color = variable)) + geom_line()
