# YASARA-log-file-analysis
## Processes YASARA log files of protein-protein docked models to extract and visualize information

### Prerequisite installations
* YASARA
* Python
* R

### Python Dependencies
* re
* pandas
* numpy
* matplotlib
* python-docx

### R Dependencies
* ggplot2
* dplyr
* reshape

### Introduction
YASARA is a molecular visualization software with many functions built into it. One such
function is the analysing interactions feature which predicts possible interactions between
different molecules, atoms, or residues. This repository was developed to better utilize
the output provided by this function for the purpose of docked-model studies. 

### Guidelines for creating log files (.docx)
* You will need to create three log files; one for contacts, one for hbonds, and one for
other_interactions. Name these files as contacts.docx, hbonds.docx, and other_interactions.docx
* Rename your post-docking pdb files using numbers. Eg: 0.pdb, 1.pdb etc...
* Include these pdb filenames as the first line of a page in the docx file.
* Copy and paste the interactions information from YASARA for each pdb file
* Insert page breaks at the end of interactions information for each pdb file

### How do I use this?
1. Make log files using YASARA and your post-docking pdb files following the guidelines provided (insert link)
2. Create a new folder inside the repository and place the three log files inside the folder
3. Run the scripts contacts_to_csv.py, hbonds_to_csv.py, and other_interactions_to_csv.py to convert log file data
into csv files (tabulates the data)
    - eg: to run the file contacts_to_csv.py, you'd enter this in the terminal:
    ```python3 contacts_to_csv.py path/to/input_folder/```
    - Each run will produce a csv file in your input folder
4. Run tables_for_line_plot.py to create contacts csv files for each pdb file.
5. Run the "line plot of contacts" R chunk from the R_analysis.Rmd notebook on each file
to generate its respective line plots
6. Run table_for_heat_diagram.py before scoring_sort.py
7. Run domain_binding_measure.py
8. Run and adjust ranges for loop_contributions.py

domain_binding_measure
loop_contributions
scoring_sort
table_for_heat_diagram
table_for_line_plot
analysis.R
