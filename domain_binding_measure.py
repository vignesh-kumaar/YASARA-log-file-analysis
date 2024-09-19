"""
Identifies the ratio of residues in the N domain vs the C domain of the binding protein.
 . specify the sequence range of the N domain using -s and -e

example command:
python3 domain_binding_measure.py -i input_files/ -s 13 -e 139
"""

import argparse
import pandas as pd
from modules import io_utils


def main():
    """ Business Logic """
    args = get_cli_args()
    infile = args.input
    start_seq = args.start
    end_seq = args.end
    h_bonds_table = pd.read_csv(infile + 'hbonds.csv', sep='\t')
    other_interactions_table = pd.read_csv(infile + 'other_interactions.csv', sep='\t')

    other_interactions_table['N.of interactions'] = pd.to_numeric(
        other_interactions_table['N.of interactions'], errors='coerce')
    other_interactions_table['N.of interactions'] = (other_interactions_table['N.of interactions']
                                                     .fillna(0))
    print(other_interactions_table)
    n_domain_count, c_domain_count = 0, 0
    for row in range(len(h_bonds_table)):
        if start_seq <= io_utils.find_sequence_number(
                h_bonds_table.at[row, 'Bacteriocin residue and number']) <= end_seq:
            n_domain_count = n_domain_count + 1
        elif h_bonds_table.at[row, 'Bacteriocin residue and number'] != '-':
            c_domain_count = c_domain_count + 1
    for row in range(len(other_interactions_table)):
        if start_seq <= io_utils.find_sequence_number(
                other_interactions_table.at[row, 'Bacteriocin residue and number']) <= end_seq:
            n_domain_count = n_domain_count + other_interactions_table.at[row, 'N.of interactions']
        else:
            c_domain_count = c_domain_count + other_interactions_table.at[row, 'N.of interactions']

    print(n_domain_count / c_domain_count)


def get_cli_args():
    """ Obtains arguments from the user """
    parser = argparse.ArgumentParser(description='Provide path+name of YASARA log file data')

    parser.add_argument('-i', '--input', type=str, help='provide the YASARA log file name',
                        required=True)
    parser.add_argument('-s', '--start', type=int, help='initial sequence number for N domain',
                        default= 13)
    parser.add_argument('-e', '--end', type=int, help='end sequence number for N domain',
                        default=139)
    return parser.parse_args()


if __name__ == "__main__":
    main()
