"""
Creates n tables for the number of contacts associated with each receptor and binding residue for the n pdb files
used in the Contacts analysis.
example terminal prompt: python3 tables_for_line_plot.py
Input: Contacts.csv
Output: line_plot_tables folder
"""


import pandas as pd
import os
import re
from modules import io_utils


def main():
    args = io_utils.get_cli_args()
    infile = args.input
    # Reads the contacts.csv file
    large_table = pd.read_csv(infile + "contacts.csv", sep="\t")

    # Uses the table as input and returns an empty dict of lists and a search range list with the
    # row indices of the table that distinguish between pdb files and the first and last index of the table.
    dict_of_list, search_range = io_utils.extract_pdb_information_from_table(large_table)

    # Writes information in the dict of lists to align with the format required to create line plots in R
    dict_of_list = write_to_dict_of_list(large_table, dict_of_list, search_range)

    # replaces receptor residue dashes with nearest previous receptor residue
    temp_receptor_residue = None
    no_dash_incomplete_list = {val: [] for val in range(len(search_range) - 1)}
    for key in range(len(search_range) - 1):
        pos = 0
        for entry in dict_of_list[key]:
            receptor_residue, binding_residue = entry[0].split('+')
            if receptor_residue == '-':
                no_dash_incomplete_list[key].append((temp_receptor_residue + '+' + binding_residue,
                                                     dict_of_list[key][pos][1]))
            else:
                temp_receptor_residue = receptor_residue
            pos = pos + 1
    # print(no_dash_incomplete_list)

    # Modification:
    final_position = find_sequence_number(no_dash_incomplete_list[len(search_range) - 2][-1][0])
    number_of_files = len(large_table[large_table['PDB file number'].str.contains('pdb', na=False)])
    pairs_contacts = {val: [] for val in range(number_of_files)}
    pdb_known_residues_final_position = []
    for key in range(len(search_range) - 1):
        pdb_known_residues_final_position.append(no_dash_incomplete_list[key][-1])
    known_final_positions = []
    for key in range(len(search_range) - 1):
        for residue in reversed(no_dash_incomplete_list[key]):
            if not residue[0].startswith('-'):
                known_final_positions.append(residue)
                break

    for key in range(len(search_range) - 1):
        null_index = 0
        repeat_index_old = None
        for position in range(len(no_dash_incomplete_list[key])):
            if position == 0:
                null_index = 1
                repeat_index = find_sequence_number(no_dash_incomplete_list[key][position][0]) - 1
                for n in range(repeat_index):
                    pairs_contacts[key].append((null_index, 0))
                    null_index = null_index + 1
                pairs_contacts[key].append(no_dash_incomplete_list[key][position])
                null_index = null_index + 1
            else:
                if repeat_index_old == 0 and repeat_index != 0:
                    null_index = null_index + 1
                repeat_index_old = repeat_index
                repeat_index = (find_sequence_number(no_dash_incomplete_list[key][position][0]) -
                                find_sequence_number(no_dash_incomplete_list[key][position - 1][0])
                                - 1)
                if repeat_index < 0:
                    repeat_index = 0
                for n in range(repeat_index):
                    # print(null_index, repeat_index_old, repeat_index)
                    pairs_contacts[key].append((null_index, 0))
                    null_index = null_index + 1
                pairs_contacts[key].append(no_dash_incomplete_list[key][position])

    for key in range(len(search_range) - 1):
        if find_sequence_number(known_final_positions[key][0]) != final_position:
            repeat_index = final_position - find_sequence_number(known_final_positions[key][0])
            for n in range(repeat_index):
                pairs_contacts[key].append((null_index, 0))
                null_index = null_index + 1

    # Uses this dict of lists to write information from each pdb in a separate csv file.
    write_output_to_csv_files(infile, pairs_contacts)


def write_to_dict_of_list(large_table, dict_of_list, search_range):
    # goes through elements associated with a pdb file per loop and appends the receptor + binding protein pair, and the
    # strength of contacts to the dict of lists
    for i, (start, end) in enumerate(zip(search_range[:-1], search_range[1:])):
        for row, column in large_table.iloc[start + 1:end].iterrows():
            dict_of_list[i].append((column.iloc[1] + '+' + column.iloc[2], column.iloc[3]))

    return dict_of_list


def write_output_to_csv_files(infile, dict_of_list):
    # creates a directory to store output
    title = infile + 'line_plot_tables/contacts_of_pdb_.csv'
    if not os.path.exists(os.path.dirname(title)):
        os.mkdir(os.path.dirname(title))

    # writes each list in dict_of_list in separate tables
    last_slash_index = title.rfind('/')
    for j in dict_of_list:
        df = pd.DataFrame(dict_of_list[j])
        specific_title = title[:last_slash_index + 17] + str(j) + title[last_slash_index + 17:]
        df.to_csv(specific_title, sep='\t', encoding='utf-8', index=False)


def find_sequence_number(residue):
    query = re.search(r'(\d+)', residue)
    return int(query.group()) if query else float('inf')


if __name__ == "__main__":
    main()
