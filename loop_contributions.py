"""


example command:
python3 loop_contributions.py -i input_files/
"""

import numpy as np
import pandas as pd
from modules import io_utils
import matplotlib.pyplot as plt


def main():
    args = io_utils.get_cli_args()
    infile = args.input

    interactions_table = pd.read_csv(infile + 'interactions_table.csv', sep='\t')
    interactions_table['Receptor residue'] = interactions_table['Receptor residue'].apply(io_utils.find_sequence_number).astype(int)

    track_particular_loops(interactions_table)


def track_particular_loops(interactions_table):
    interacting_residues_count = len(interactions_table[interactions_table['Contacts across PDBs'] >= 4])
    interactions_table = interactions_table[interactions_table['Contacts across PDBs'] >= 4]
    print(interactions_table)

    # L 4,6,7 for pseudomonas fluorescens bamA
    # range_1 = other_interactions_table['Receptor residue'].between(523, 538)
    # range_2 = other_interactions_table['Receptor residue'].between(556, 570)
    # range_3 = other_interactions_table['Receptor residue'].between(616, 630)
    # range_4 = other_interactions_table['Receptor residue'].between(691, 702)
    # range_5 = other_interactions_table['Receptor residue'].between(712, 725)
    # range_6 = other_interactions_table['Receptor residue'].between(753, 765)
    # condition = (range_1) | (range_2) | (range_3) | (range_4) | (range_5) | (range_6)

    # L 4, 6, 7 including ELs for pseudomonas fluorescens bamA
    # range_1 = interactions_table['Receptor residue'].between(523, 570)
    # range_2 = interactions_table['Receptor residue'].between(616, 702)
    # range_3 = interactions_table['Receptor residue'].between(712, 765)
    # condition = range_1 | range_2 | range_3

    # L 4, 6, 7 including ELs for ecoli bamA
    range_1 = interactions_table['Receptor residue'].between(522, 579)
    range_2 = interactions_table['Receptor residue'].between(627, 720)
    range_3 = interactions_table['Receptor residue'].between(734, 778)
    condition = range_1 | range_2 | range_3

    df = interactions_table[condition]
    df = df[df['Contacts across PDBs'] > 0]
    print(df)
    particular_loop_count = len(df)
    print(particular_loop_count, interacting_residues_count)
    interactions_table.plot(x='Receptor residue',  y='Contacts across PDBs', style=['o'])
    plt.show()


if __name__ == "__main__":
    main()
