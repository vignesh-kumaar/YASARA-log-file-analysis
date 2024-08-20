pip install python-docx # for reading doc files
pdb file names should be named in numbers and entered as such in the log files. Eg: 0.pdb, 1.pdb, etc...

steps to use:
1. paste log files in input files folder or modify path accordingly
2. Use contacts_to_csv.py, hbonds_to_csv.py, and other_interactions_to_csv.py to convert log file data
into csv files
3. Use table_for_heat_diagram.py, and tables_for_line_plot.py to create a heat map from the contacts.csv file