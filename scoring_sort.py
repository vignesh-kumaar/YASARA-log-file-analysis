"""
python3 scoring_sort.py -i input_files/
python3 scoring_sort.py -i ecoli_bamA_mosselii_llpA1/analysis/
"""

import pandas as pd
import re
import argparse


def main():
    args = get_cli_args()
    infile = args.input

    # Extracts data from csv files
    # contacts_table = pd.read_csv("contacts.csv", sep='\t')
    h_bonds_table = pd.read_csv(infile + "hbonds.csv", sep='\t')
    other_interactions_table = pd.read_csv(infile + "other_interactions.csv", sep='\t')
    contacts_unsorted_table = pd.read_csv(infile + "number_of_contacts_with_PDBs.csv", sep='\t')
    interactions_table = contacts_unsorted_table
    other_interactions_table['Interaction strength'] = pd.to_numeric(other_interactions_table['Interaction strength'],
                                                                     errors="coerce")
    print(other_interactions_table)

    # write the interactions into the interactions_table
    interactions_table['h_bonds'] = ''

    # Add information into interactions_table
    for i in range(1, len(contacts_unsorted_table) + 1):
        res_count = 0
        # using h_bonds_table:
        for j in range(len(h_bonds_table)):
            if 'pdb' in h_bonds_table.at[j, 'PDB file number']:
                continue
            # elif h_bonds_table.at[j, 'Receptor residue and number'] != '-':
            elif i == find_sequence_number(h_bonds_table.at[j, 'Receptor residue and number']):
                res_count = res_count + 1
                while h_bonds_table.at[j, 'Receptor residue and number'] == '-':
                    if 'pdb' in h_bonds_table.at[j, 'PDB file number']:
                        break
                    res_count = res_count + 1
                    j = j + 1
            else:
                continue
        interactions_table.at[i-1, 'h_bonds'] = res_count

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
    interactions_table[['Ionic', 'Hydrophobic', 'CationPi', 'PiPi']] = (
        interactions_table[['Ionic', 'Hydrophobic', 'CationPi', 'PiPi']].fillna(0))
    interactions_table = interactions_table.sort_values(by='Contacts across PDBs', ascending=False)
    print(interactions_table)
    interactions_table.to_csv(infile + 'interactions_table.csv', sep='\t', encoding='utf-8', index=False)


def find_sequence_number(residue):
    query = re.search(r'\d+$', residue)
    return int(query.group()) if query else float('inf')


def get_cli_args():
    parser = argparse.ArgumentParser(description='Provide path to YASARA log file data')

    parser.add_argument('-i', '--input', type=str, help='provide the path to YASARA log files', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    main()
