"""
File: YASARA_data_to_csv.py

This program takes YASARA log file data and outputs necessary information into a csv file

example command: python3 hbonds_to_csv.py -i input_files/
python3 hbonds_to_csv.py -i ecoli_bamA_mosselii_llpA1/analysis/
"""

from modules import io_utils
from docx import Document
import pandas as pd


def main():
    """Business Logic"""

    # declaring and initializing data
    args = io_utils.get_cli_args()
    infile = args.input

    # reading log data into python
    df = read_hbonds_to_df(infile + 'hbonds.docx')

    # writing to output file
    df.to_csv(infile + 'hbonds.csv', sep='\t', encoding='utf-8', index=False)


def read_hbonds_to_df(file=None):
    # Step 1
    # Initialize empty list and open log file using Document method
    data_array = []
    logs = Document(file)

    # Iterate through every line of the log file and extract required data into the list
    for line in logs.paragraphs:
        if line.text:
            if len(line.text) >= 3:
                if line.text[2] == "p" or line.text[3] == "p":
                    data_array.append(line.text)
            words = line.text.split()
            if words[0] == "Residue":
                data_array.append(words[1] + words[3])
            elif words[0] == "Atom":
                data_array.append(words[1:2] + [words[2] + words[3]] + words[5:6] + words[9:10] + [words[10] + words[11]] + words[17:18] + words[22:23])
            # elif len(words) >= 2:
            #    if words[1] == "hydrogen":
            #        data_array.append(words[16])

    # Convert list into df for convenience in the next step
    df = pd.DataFrame(data_array)
    df = _format_df(df)
    return df


def _format_df(df):
    # Step 2
    # Initialize variable and new df for the next step
    receptor_residue_to_add = None
    output_df = pd.DataFrame(columns=[
        'PDB file number',
        'Receptor residue and number',
        'Bacteriocin residue and number',
        'Receptor atom/group',
        'Accepts/Rejects',
        'Bacteriocin atom/group',
        'Bond length',
        'Bond energy'
    ])

    # Iterate through rows of df and store bacteriocin contact residues as a separate index associated
    # with a given receptor contact residue
    for index, row in df.iterrows():
        value = row[0]
        if isinstance(value, str):  # Check if the entry is a non-list (PDB/receptor residue)
            if "pdb" in value:
                new_row = pd.DataFrame([{
                    'PDB file number': value,
                    'Receptor residue and number': '-',
                    'Bacteriocin residue and number': '-',
                    'Receptor atom/group': '-',
                    'Accepts/Rejects': '-',
                    'Bacteriocin atom/group': '-',
                    'Bond length': '-',
                    'Bond energy': '-'
                }])
                output_df = pd.concat([output_df, new_row], ignore_index=True)
            else:
                receptor_residue_to_add = value
        elif isinstance(value, list):
            new_row = pd.DataFrame([{
                'PDB file number': '-',
                'Receptor residue and number': receptor_residue_to_add,
                'Bacteriocin residue and number': value[4],
                'Receptor atom/group': value[0],
                'Accepts/Rejects': value[2],
                'Bacteriocin atom/group': value[3],
                'Bond length': value[5],
                'Bond energy': value[6]
            }])
            # Concatenate the new row DataFrame to the output DataFrame
            output_df = pd.concat([output_df, new_row], ignore_index=True)
            # After the first addition, set receptor_residue_to_add to None for subsequent additions
            receptor_residue_to_add = None

    # Replace None with '-' and return output df
    output_df.fillna('-', inplace=True)
    print(output_df)
    return output_df


if __name__ == "__main__":
    main()
