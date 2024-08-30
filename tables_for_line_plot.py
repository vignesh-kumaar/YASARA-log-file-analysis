"""
Creates n tables with contacts information for the n pdb files used in the YASARA analysis

example command:
python3 tables_for_line_plot.py -i input_files/

Input: contacts.csv
Output: line_plot_tables folder
"""

import os
import re
import pandas as pd

from modules import io_utils


def main():
    """ Business Logic """

    args = io_utils.get_cli_args()
    infile = args.input
    # Reads the contacts.csv file
    contacts_table = pd.read_csv(infile + "contacts.csv", sep="\t")

    # Uses the table as input and returns an empty dict of lists and a search range list with the
    # row indices of the table that separate pdb files and the first and last index of the table.
    dict_of_contacts, search_range = io_utils.extract_pdb_information_from_table(contacts_table)

    # Formats dict_of_contacts to be suitable for creation of line plots in R
    dict_of_contacts = write_to_dict_of_contacts(contacts_table, dict_of_contacts, search_range)

    # replaces receptor residue dashes with nearest previous receptor residue
    no_dash_incomplete_dict = fill_dashes_with_receptors(search_range, dict_of_contacts)

    # writes interacting contact pairs to csv files
    write_contact_pairs_to_csv(no_dash_incomplete_dict, infile)
    
    # To fill the distance between adjacent interacting residues with zeros:
    contacts_grid = fill_zeros_for_non_interactions(no_dash_incomplete_dict, search_range,
                                                     contacts_table)

    # writes contact grid information to csv files
    write_output_to_csv_files(infile, contacts_grid)


def write_to_dict_of_contacts(contacts_table, dict_of_contacts, search_range):
    """
    Adds information to a dict_of_contacts suitable to create line plots with.
    :param contacts_table: data from the contacts.csv file generated using contacts_to_csv.py
    :param dict_of_contacts: empty dict of lists to be filled
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :return: the dictionary of lists dict_of_contacts
    """

    # Goes through elements associated with a pdb file per loop and appends the receptor
    # + binding protein pair, and the strength of contacts to dict_of_contacts
    for i, (start, end) in enumerate(zip(search_range[:-1], search_range[1:])):
        for _, column in contacts_table.iloc[start + 1:end].iterrows():
            dict_of_contacts[i].append((column.iloc[1] + '+' + column.iloc[2], column.iloc[3]))

    return dict_of_contacts


def fill_dashes_with_receptors(search_range, dict_of_contacts):
    """
    dict_of_contacts currently contains '-' as entries when previous receptor residue sequences
    are repeated. For the purpose of line_plots,  all entries are required to include the receptor
    residue sequence. The function generates a new variable that fits the above needs.
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :param dict_of_contacts: dictionary of lists with '-' entries that need to be modified
    :return: dictionary of lists no_dash_incomplete_dict
    """

    temp_receptor_residue = None
    no_dash_incomplete_dict = {val: [] for val in range(len(search_range) - 1)}
    for key in range(len(search_range) - 1):
        pos = 0
        for entry in dict_of_contacts[key]:
            receptor_residue, binding_residue = entry[0].split('+')
            if receptor_residue == '-':
                no_dash_incomplete_dict[key].append((temp_receptor_residue + '+' + binding_residue,
                                                     dict_of_contacts[key][pos][1]))
            else:
                temp_receptor_residue = receptor_residue
            pos = pos + 1
    return no_dash_incomplete_dict


def write_contact_pairs_to_csv(no_dash_incomplete_dict, infile):
    # creates a directory to store output
    title = infile + 'non_zero_contacts_tables/contacts_of_pdb_.csv'
    if not os.path.exists(os.path.dirname(title)):
        os.mkdir(os.path.dirname(title))

    # writes each list in no_dash_incomplete_dict in separate tables
    last_slash_index = title.rfind('/')
    for j in no_dash_incomplete_dict:
        df = pd.DataFrame(no_dash_incomplete_dict[j])
        specific_title = title[:last_slash_index + 17] + str(j) + title[last_slash_index + 17:]
        df.to_csv(specific_title, sep='\t', encoding='utf-8', index=False)


def fill_zeros_for_non_interactions(no_dash_incomplete_dict, search_range, contacts_table):
    """
    To show an accurate representation of distances between interacting residues, null entries
    need to be filled in between them.
    :param no_dash_incomplete_dict: dict of lists with no distances between interacting residues
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :param contacts_table: data from the contacts.csv file generated using contacts_to_csv.py
    :return: dictionary of lists contacts_grid
    """

    # Declaring and initializing variables
    number_of_files = len(contacts_table[contacts_table['PDB file number'].str.contains(
        'pdb', na=False)])
    contacts_grid = {val: [] for val in range(number_of_files)}
    pdb_known_residues_final_position = []

    # Record positions of final interacting residues in each pdb file and store the position
    # for the final pdb file to use later
    for key in range(len(search_range) - 1):
        pdb_known_residues_final_position.append(no_dash_incomplete_dict[key][-1])
    final_position = find_sequence_number(pdb_known_residues_final_position[-1][0])

    # Fill zeros to account for distance between adjacent interacting residues
    for key in range(len(search_range) - 1):
        null_index = 0
        repeat_index_old = None
        for position in range(len(no_dash_incomplete_dict[key])):
            if position == 0:
                null_index = 1
                repeat_index = find_sequence_number(no_dash_incomplete_dict
                                                             [key][position][0]) - 1
                for _ in range(repeat_index):
                    contacts_grid[key].append((null_index, 0))
                    null_index = null_index + 1
                contacts_grid[key].append(no_dash_incomplete_dict[key][position])
                null_index = null_index + 1
            else:
                if repeat_index_old == 0 and repeat_index != 0:
                    null_index = null_index + 1
                repeat_index_old = repeat_index
                repeat_index = (find_sequence_number(no_dash_incomplete_dict
                                                              [key][position][0]) -
                                find_sequence_number(no_dash_incomplete_dict
                                                              [key][position - 1][0]) - 1)
                repeat_index = max(repeat_index, 0)
                for _ in range(repeat_index):
                    # print(null_index, repeat_index_old, repeat_index)
                    contacts_grid[key].append((null_index, 0))
                    null_index = null_index + 1
                contacts_grid[key].append(no_dash_incomplete_dict[key][position])
    # Fill zeros between the final interacting residue and the terminal end of the protein.
    for key in range(len(search_range) - 1):
        if find_sequence_number(pdb_known_residues_final_position[key][0]) != final_position:
            repeat_index = final_position - find_sequence_number(
                pdb_known_residues_final_position[key][0])
            for _ in range(repeat_index):
                contacts_grid[key].append((null_index, 0))
                null_index = null_index + 1
    return contacts_grid

def write_output_to_csv_files(infile, contacts_grid):
    """
    Writes data of each PDB file from contacts_grid in different csv files
    :param infile: path to the I/O directory provided by the user
    :param contacts_grid: dict of list containing formatted contacts information
    :return: None
    """

    # creates a directory to store output
    title = infile + 'line_plot_tables/contacts_of_pdb_.csv'
    if not os.path.exists(os.path.dirname(title)):
        os.mkdir(os.path.dirname(title))

    # writes each list in dict_of_contacts in separate tables
    last_slash_index = title.rfind('/')
    for j in contacts_grid:
        df = pd.DataFrame(contacts_grid[j])
        specific_title = title[:last_slash_index + 17] + str(j) + title[last_slash_index + 17:]
        df.to_csv(specific_title, sep='\t', encoding='utf-8', index=False)


def find_sequence_number(entry):
    """
    Given a string, returns the first instance of a number before a plus sign. Practically, it
    captures the sequence number of the receptor residue in a given 'receptor residue + binding
    residue' interaction entry.
    :param entry: string with a receptor residue and sometimes a binding residue separated by '+'.
    :return: integer containing the receptor residue sequence number or 'inf'
    """

    query = re.search(r'(\d+)', entry)
    return int(query.group()) if query else float('inf')


if __name__ == "__main__":
    main()
