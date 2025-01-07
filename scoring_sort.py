"""
Sorts receptor residues as likely binding sites using information from all identified interactions.

NOTE: Run table_for_heat_diagram.py before running this script
example command:
python3 scoring_sort.py -i input_files/
"""

import pandas as pd
import re
import argparse


def main():
    """ Business Logic """
    args = get_cli_args_scoring_sort()
    infile = args.input1
    sorting_type = args.input2

    interactions_table = fill_interactions_table(infile)

    # Replace blank/none entries with 0s
    interactions_table[['Ionic', 'Hydrophobic', 'CationPi', 'PiPi']] = (
        interactions_table[['Ionic', 'Hydrophobic', 'CationPi', 'PiPi']].fillna(0))

    if sorting_type == 1:
        # To display interactions in order of sequence number
        interactions_table['sequence_number'] = interactions_table['Receptor residue'].apply(find_sequence_number)
        interactions_table = interactions_table.sort_values(by='sequence_number', ascending=True)
        interactions_table = interactions_table.drop(columns=['sequence_number'])
    elif sorting_type == 2:
        # To display interactions based on hydrogen bonding
        interactions_table = interactions_table.sort_values(by=['h_bonds', 'mean_h_bond_energy'], ascending=False)
    elif sorting_type == 3:
        interactions_table['mean_h_bond_energy'] = interactions_table['mean_h_bond_energy'].astype(float)
        interactions_table = interactions_table.sort_values(by=['mean_h_bond_energy'], ascending=False)

    print(interactions_table)

    # Write interactions_table to output
    interactions_table.to_csv(infile + 'interactions_table.csv', sep='\t', encoding='utf-8', index=False)

    # Create h_bonds_matrix and write to output
    h_bonds_matrix = interactions_table[['Receptor residue', 'h_bonds', 'mean_h_bond_energy']]
    h_bonds_matrix = h_bonds_matrix[h_bonds_matrix['h_bonds'] > 0]
    print(h_bonds_matrix)
    h_bonds_matrix.to_csv(infile + 'h_bonds_matrix.csv', sep='\t', encoding='utf-8', index=False)


def fill_interactions_table(infile):
    # Extracts data from csv files
    h_bonds_table = pd.read_csv(infile + "hbonds.csv", sep='\t')
    other_interactions_table = pd.read_csv(infile + "other_interactions.csv", sep='\t')
    contacts_unsorted_table = pd.read_csv(infile + "number_of_contacts_with_PDBs.csv", sep='\t')
    interactions_table = contacts_unsorted_table
    other_interactions_table['Interaction strength'] = pd.to_numeric(other_interactions_table['Interaction strength'],
                                                                     errors="coerce")

    # Create new columns in the interactions_table
    interactions_table['h_bonds'] = ''
    # interactions_table['h_bond_energies'] = ''
    interactions_table['mean_h_bond_energy'] = ''

    # Add information into interactions_table
    for i in range(1, len(contacts_unsorted_table) + 1):
        h_bond_energies = []
        res_count = 0
        # using h_bonds_table:
        for j in range(len(h_bonds_table)):
            if 'pdb' in h_bonds_table.at[j, 'PDB file number']:
                continue
            elif i == find_sequence_number(h_bonds_table.at[j, 'Receptor residue and number']):
                res_count = res_count + 1
                h_bond_energies.append(float(h_bonds_table.at[j, 'Bond energy']))
                while h_bonds_table.at[j, 'Receptor residue and number'] == '-':
                    if 'pdb' in h_bonds_table.at[j, 'PDB file number']:
                        break
                    res_count = res_count + 1
                    j = j + 1
            else:
                continue
        interactions_table.at[i-1, 'h_bonds'] = res_count
        # interactions_table.at[i-1, 'h_bond_energies'] = h_bond_energies
        if interactions_table.at[i-1, 'h_bonds'] > 0:
            interactions_table.at[i-1, 'mean_h_bond_energy'] = (sum(h_bond_energies) /
                                                                len(h_bond_energies))
        else:
            interactions_table.at[i-1, 'mean_h_bond_energy'] = 0

        # using other_interactions_table:
        (ionic_interaction_strength, hydrophobic_interaction_strength,
         cation_pi_interaction_strength, pi_pi_interaction_strength) = 0, 0, 0, 0
        for k in range(len(other_interactions_table)):
            if 'pdb' in other_interactions_table.at[k, 'PDB file number']:
                continue
            elif i == find_sequence_number(other_interactions_table.at[k, 'Receptor residue and number']):
                if other_interactions_table.at[k, 'Type of Interaction'] == 'Ionic':
                    ionic_interaction_strength = (ionic_interaction_strength +
                                                  other_interactions_table.at[k, 'Interaction strength'])
                elif other_interactions_table.at[k, 'Type of Interaction'] == 'Hydrophobic':
                    hydrophobic_interaction_strength = (hydrophobic_interaction_strength +
                                                        other_interactions_table.at[k, 'Interaction strength'])
                elif other_interactions_table.at[k, 'Type of Interaction'] == 'CationPi':
                    cation_pi_interaction_strength = (cation_pi_interaction_strength +
                                                      other_interactions_table.at[k, 'Interaction strength'])
                elif other_interactions_table.at[k, 'Type of Interaction'] == 'PiPi':
                    pi_pi_interaction_strength = (pi_pi_interaction_strength +
                                                  other_interactions_table.at[k, 'Interaction strength'])
                else:
                    continue
                while other_interactions_table.at[k, 'Receptor residue and number'] == '-':
                    if 'pdb' in other_interactions_table.at[k, 'PDB file number']:
                        break

                    if other_interactions_table.at[k, 'Type of Interaction'] == 'Ionic':
                        ionic_interaction_strength = (ionic_interaction_strength +
                                                      other_interactions_table.at[k, 'Interaction strength'])
                    elif other_interactions_table.at[k, 'Type of Interaction'] == 'Hydrophobic':
                        hydrophobic_interaction_strength = (hydrophobic_interaction_strength +
                                                            other_interactions_table.at[k, 'Interaction strength'])
                    elif other_interactions_table.at[k, 'Type of Interaction'] == 'CationPi':
                        cation_pi_interaction_strength = (cation_pi_interaction_strength +
                                                          other_interactions_table.at[k, 'Interaction strength'])
                    elif other_interactions_table.at[k, 'Type of Interaction'] == 'PiPi':
                        pi_pi_interaction_strength = (pi_pi_interaction_strength +
                                                      other_interactions_table.at[k, 'Interaction strength'])
                    else:
                        continue
                    k = k + 1
            else:
                continue
            interactions_table.at[i - 1, 'Ionic'] = ionic_interaction_strength
            interactions_table.at[i - 1, 'Hydrophobic'] = hydrophobic_interaction_strength
            interactions_table.at[i - 1, 'CationPi'] = cation_pi_interaction_strength
            interactions_table.at[i - 1, 'PiPi'] = pi_pi_interaction_strength

    return interactions_table


def find_sequence_number(residue):
    query = re.search(r'\d+$', residue)
    return int(query.group()) if query else float('inf')


def get_cli_args_scoring_sort():
    parser = argparse.ArgumentParser(description='Provide path to YASARA log file data and sorting type')
    parser.add_argument('-i1', '--input1', type=str,
                        help='provide the path to YASARA log files', required=True)
    parser.add_argument('-i2','--input2', type=int,
                        help='choose between 1."sequence number", 2."file occurences + h_bond_energy"'
                                                         'and 3."h_bond_energy"')
    args = parser.parse_args()
    # args.input1 = input(f"Provide path to YASARA log file data and sorting type")
    # args.input2 = input(f"choose between 1.'sequence number', 2.'file occurences + h_bond_energy'"
    #                    f" and 3.'h_bond_energy'")
    if args.input1 is None:
        args.input1 = input(f"Please provide a path to the log file")
    if args.input2 is None:
        args.input2 = int(input(f"Please provide a sorting type between\n1. 'sequence number'\t"
              f"2.file occurences + h_bond_energy\t3.h_bond_energy\n"))
    return args


if __name__ == "__main__":
    main()
