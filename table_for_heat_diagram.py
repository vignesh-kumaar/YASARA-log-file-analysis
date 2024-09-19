"""
Generates a heat diagram with receptor residue vs PDB file number for the axes and colored
based on the strength of the interaction.

example command:
python3 table_for_heat_diagram.py -i input_files/
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from modules import io_utils


def main():
    """ Business Logic """

    args = io_utils.get_cli_args()
    infile = args.input

    # Reads the contacts.csv file
    contacts_table = pd.read_csv(infile + "contacts.csv", sep="\t")

    # extracts an empty dict of lists and a search range with the rows corresponding
    # to a change in pdb files in the table
    pdb_dict_residue_pairs, search_range = (io_utils.extract_pdb_information_from_table
                                            (contacts_table))

    # fills the dict of lists with residue pair information of pairs with non-zero contacts
    pdb_dict_residue_pairs = fill_residues_non_zero_contacts(contacts_table,
                                                             pdb_dict_residue_pairs, search_range)

    # To fill null entries for row names:
    complete_residue_list, final_position = fill_row_names(pdb_dict_residue_pairs)

    # Filling known contacts:
    contacts_matrix = fill_known_contacts(contacts_table, search_range)

    # Finding sorted receptor residues with non-zero contacts from the residue pairs:
    identified_receptor_residues_individual_list = (fill_sorted_receptor_residues
                                                    (pdb_dict_residue_pairs, search_range))

    # Finding distances between each adjacent residues in every pdb file and filling
    # in 0's for a complete grid:
    complete_contacts_matrix = fill_contacts_matrix(contacts_matrix,
                                                    identified_receptor_residues_individual_list,
                                                    search_range, final_position)

    # side project: sort residues by number of contacts
    sort_receptor_residues_by_contacts(complete_contacts_matrix, complete_residue_list, infile)

    # Formats the contacts matrix and residues list
    df, complete_residue_list = format_matrix_and_row_labels(complete_contacts_matrix,
                                                             complete_residue_list)

    # Create Heat Map:
    create_heat_map(df, complete_residue_list, search_range, infile)


def fill_residues_non_zero_contacts(contacts_table, pdb_dict_residue_pairs, search_range):
    """
    Fills empty dict of lists with all contact pairs
    :param contacts_table: data from the contacts.csv file generated using contacts_to_csv.py
    :param pdb_dict_residue_pairs: empty dict of lists to be filled with contact pair information
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :return: dictionary of lists pdb_dict_residue_pairs
    """

    for i, (start, end) in enumerate(zip(search_range[:-1], search_range[1:])):
        for _, column in contacts_table.iloc[start+1:end].iterrows():
            pdb_dict_residue_pairs[i].append(column.iloc[1] + '+' + column.iloc[2])
    return pdb_dict_residue_pairs


def fill_row_names(pdb_dict_residue_pairs):
    """
    Creates a complete list of residues for the receptor protein. All interacting residues are
    stored with the three letter amino acid code and sequence positions whereas non-interacting
    residues are stored with only their sequence number.
    :param pdb_dict_residue_pairs: dict of list containing contact pair information
    :return: list containing labels for each residue of the receptor and an integer
    containing the terminal end of the receptor protein.
    """

    identified_receptor_residues_set = set()
    for key in pdb_dict_residue_pairs:
        for element in pdb_dict_residue_pairs[key]:
            residue = element.split('+')[0]
            if '-' not in residue:
                identified_receptor_residues_set.add(residue)

    identified_receptor_residues_list = list(identified_receptor_residues_set)
    sorted_receptor_residues = sorted(identified_receptor_residues_list,
                                      key=io_utils.find_sequence_number)
    # sorted_receptor_residues now contains the non-zero contact row headings for our heat map
    identified_sequence_numbers = [int(element[3:]) for element in sorted_receptor_residues]
    final_position = identified_sequence_numbers[-1]
    tally = 0
    complete_residue_list = []
    for i, (start, end) in enumerate(zip(identified_sequence_numbers[:-1],
                                         identified_sequence_numbers[1:])):
        if i == 0:
            [complete_residue_list.append(str(n)) for n in range(1, start)]
        while identified_sequence_numbers[tally] <= start:
            complete_residue_list.append(sorted_receptor_residues[tally])
            tally = tally + 1
        [complete_residue_list.append(str(n)) for n in range(start+1, end)]
        if end == final_position:
            complete_residue_list.append(sorted_receptor_residues[tally])
    return complete_residue_list, final_position


def fill_known_contacts(contacts_table, search_range):
    """
    Stores contact values in a list of list
    :param contacts_table: data from the contacts.csv file generated using contacts_to_csv.py
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :return: dictionary of lists pdb_dict_residue_pairs
    :return: list of lists contacts_matrix
    """

    contacts_matrix = [[] for n in range(len(search_range) - 1)]
    pdb_number = -1
    for row in range(len(contacts_table)):
        if 'pdb' in contacts_table['PDB file number'][row]:
            pdb_number = pdb_number + 1
        elif contacts_table['Receptor residue and number'][row] != '-':
            contacts_matrix[pdb_number].append(_sum_of_contacts_for_residue(contacts_table, row))
        else:
            continue
    return contacts_matrix


def _sum_of_contacts_for_residue(contacts_table, row):
    """
    sums contact values for a given receptor residue
    :param contacts_table: data from the contacts.csv file generated using contacts_to_csv.py
    :param row: row number with the first entry of the given receptor residue in contacts_table
    :return: integer with the sum of contacts for the given receptor residue
    """

    contacts_count = int(contacts_table['Strength of contacts'][row])
    row = row + 1
    if row >= len(contacts_table):
        return contacts_count
    while contacts_table['Receptor residue and number'][row] == '-':
        if 'pdb' in contacts_table['PDB file number'][row]:
            break
        contacts_count = contacts_count + int(contacts_table['Strength of contacts'][row])
        row = row + 1
        if row >= len(contacts_table):
            break
    return contacts_count


def fill_sorted_receptor_residues(pdb_dict_residue_pairs, search_range):
    """
    creates a sorted list of all the interacting receptor residues from the residue pairs
    :param pdb_dict_residue_pairs: dict of list containing contact pair information
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :return: list containing all interacting receptor residues in ascending order
    """

    identified_receptor_residues_individual = [set() for n in range(len(search_range) - 1)]
    for key in pdb_dict_residue_pairs:
        for element in pdb_dict_residue_pairs[key]:
            residue = element.split('+')[0]
            if '-' not in residue:
                identified_receptor_residues_individual[key].add(residue)

    identified_receptor_residues_individual_list = [[] for n in range(len(search_range) - 1)]
    for number in range(len(search_range) - 1):
        identified_receptor_residues_individual_list[number] = sorted(list(
            identified_receptor_residues_individual[number]), key=io_utils.find_sequence_number)
    return identified_receptor_residues_individual_list


def fill_contacts_matrix(contacts_matrix, identified_receptor_residues_individual_list,
                         search_range, final_position):
    """
    The current contacts_matrix contains all non-zero contacts information across PDB files. This
    function fills zeros between two interacting receptor residues to show the distance between
    these regions in the heat plot.
    :param contacts_matrix: list of list containing non-zero contact values
    :param identified_receptor_residues_individual_list: list containing all interacting receptor
    residues in ascending order
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :param final_position: integer containing the final residue position of the receptor
    :return: List of list complete_contacts_matrix
    """

    complete_contacts_matrix = [[] for n in range(len(search_range) - 1)]
    pdb_known_residues_final_position = []
    for key in range(len(search_range) - 1):
        pdb_known_residues_final_position.append(
            identified_receptor_residues_individual_list[key][-1])

    for key in range(len(search_range) - 1):
        for position in range(len(identified_receptor_residues_individual_list[key])):
            if position == 0:
                repeat_index = io_utils.find_sequence_number(
                    identified_receptor_residues_individual_list[key][position]) - 1
                for n in range(repeat_index):
                    complete_contacts_matrix[key].append(0)
                complete_contacts_matrix[key].append(contacts_matrix[key][position])
            else:
                repeat_index = (io_utils.find_sequence_number(
                    identified_receptor_residues_individual_list[key][position])
                                - io_utils.find_sequence_number(
                            identified_receptor_residues_individual_list[key][position - 1])
                                - 1)
                for n in range(repeat_index):
                    complete_contacts_matrix[key].append(0)
                complete_contacts_matrix[key].append(contacts_matrix[key][position])
    for key in range(len(search_range) - 1):
        if io_utils.find_sequence_number(pdb_known_residues_final_position[key]) != final_position:
            repeat_index = final_position - io_utils.find_sequence_number(
                pdb_known_residues_final_position[key])
            for n in range(repeat_index):
                complete_contacts_matrix[key].append(0)
    return complete_contacts_matrix


def sort_receptor_residues_by_contacts(complete_contacts_matrix, complete_residue_list, infile):
    """
    Sorts in descending order all receptor residues based on the number of contacts with the
    binding protein.
    :param complete_contacts_matrix: list of list containing all contacts values including zeros.
    Each sublist of the list corresponds to a particular PDB file.
    :param complete_residue_list: list of all receptor residues.
    :param infile: path to the I/O directory provided by the user
    :return: None
    """

    df = pd.DataFrame(complete_contacts_matrix)
    df = df.transpose()
    row_headings = pd.Series(complete_residue_list)
    non_zero_contacts_df = df.astype(bool).sum(axis=1)
    non_zero_contacts_df = non_zero_contacts_df.to_frame()
    non_zero_contacts_df['Receptor residue'] = row_headings.values
    cols = non_zero_contacts_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    non_zero_contacts_df = non_zero_contacts_df[cols]
    non_zero_contacts_df.columns = ['Receptor residue', 'Contacts across PDBs']
    sorted_non_zero_contacts_df = non_zero_contacts_df.sort_values(by='Contacts across PDBs',
                                                                   ascending=False)
    sorted_non_zero_contacts_df.to_csv(infile + 'sorted_number_of_contacts_with_PDBs.csv', sep='\t',
                                       encoding='utf-8', index=False)
    non_zero_contacts_df.to_csv(infile + 'number_of_contacts_with_PDBs.csv', sep='\t',
                                encoding='utf-8', index=False)


def format_matrix_and_row_labels(complete_contacts_matrix, complete_residue_list):
    """
    Formats the contacts matrix and reverses the residue list. The new contacts matrix is a list of
    list with each sublist corresponding to a given residue instead of a given PDB file.
    :param complete_contacts_matrix: list of list containing all contacts values including zeros.
    Each sublist of the list corresponds to a particular PDB file.
    :param complete_residue_list: List of all receptor residues
    :return: dataframe df containing the contacts matrix and list complete_residue_list containing
    all receptor residues in descending order.
    """

    complete_contacts_matrix = list(map(list, zip(*complete_contacts_matrix)))
    df = pd.DataFrame(complete_contacts_matrix)
    # df2 = df.iloc[700:792]
    complete_residue_list = complete_residue_list[::-1]
    return df, complete_residue_list


def create_heat_map(df, complete_residue_list, search_range, infile):
    """
    Creates the heat map.
    :param df: Dataframe containing the contacts values of all receptor residue interactions
    :param complete_residue_list: List of all receptor residues in descending order
    :param search_range: list with row indices of contacts_table where a new pdb file is recorded
    :param infile: path to the I/O directory provided by the user
    :return: None
    """

    # plt.yticks(range(0, len(complete_residue_list)), complete_residue_list, fontsize=3)
    y_labels = np.arange(0, len(complete_residue_list), 10)
    x_labels = [f"{i}.pdb" for i in range(len(search_range))]
    plt.yticks(y_labels, fontsize=3)
    plt.xticks(range(len(search_range)), labels=x_labels, fontsize=4)
    plt.imshow(df, cmap='magma', interpolation='nearest', aspect=0.025)
    # aspect='auto'
    plt.gca().invert_yaxis()
    plt.tight_layout()
    df.to_csv(infile + "contacts_heat_matrix.csv", sep='\t', encoding='utf-8', index=False)
    plt.savefig(infile + "heatmap.png", dpi=1000)
    # plt.show()


if __name__ == "__main__":
    main()
