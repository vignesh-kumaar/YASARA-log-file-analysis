"""
python3 domain_binding_measure.py -i input_files/
"""

import pandas as pd
from modules import io_utils


def main():
    args = io_utils.get_cli_args()
    infile = args.input
    h_bonds_table = pd.read_csv(infile + 'hbonds.csv', sep='\t')
    other_interactions_table = pd.read_csv(infile + 'other_interactions.csv', sep='\t')

    # print(h_bonds_table)
    # print(other_interactions_table)
    other_interactions_table['N.of interactions'] = pd.to_numeric(other_interactions_table['N.of interactions'],
                                                                  errors='coerce')
    other_interactions_table['N.of interactions'] = other_interactions_table['N.of interactions'].fillna(0)
    n_domain_count, c_domain_count = 0, 0
    for row in range(len(other_interactions_table)):
        if 13 <= io_utils.find_sequence_number(
                other_interactions_table.at[row, 'Bacteriocin residue and number']) <= 139:
            n_domain_count = n_domain_count + other_interactions_table.at[row, 'N.of interactions']
        else:
            c_domain_count = c_domain_count + other_interactions_table.at[row, 'N.of interactions']

    print(n_domain_count/c_domain_count)


if __name__ == "__main__":
    main()
