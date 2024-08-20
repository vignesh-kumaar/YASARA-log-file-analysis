"""
File: YASARA_data_to_csv.py

This program takes YASARA log file data and outputs necessary information into a csv file

example command: python3 other_interactions_to_csv.py -i input_files/
python3 other_interactions_to_csv.py -i ecoli_bamA_mosselii_llpA1/analysis/
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
    df = read_other_interactions_to_df(infile + 'other_interactions.docx')

    # writing to output file
    df.to_csv(infile + 'other_interactions.csv', sep='\t', encoding='utf-8', index=False)


def read_other_interactions_to_df(file=None):
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
            if words[0] in ("Hydrophobic", "PiPi", "CationPi", "Ionic"):
                data_array.append(words[0])
            elif words[0] == "Residue":
                data_array.append(words[1] + words[3])
            elif words[0] == "to":
                data_array.append([words[2] + words[4]] + words[6:7] + words[10:11])

    # Convert list into df for convenience in the next step
    df = pd.DataFrame(data_array)
    df = _format_df(df)
    return df


def _format_df(df):
    # Step 2
    # Initialize variable and new df for the next step
    receptor_residue_to_add = None
    type_of_interaction = None
    output_df = pd.DataFrame(columns=[
        'PDB file number',
        'Receptor residue and number',
        'Bacteriocin residue and number',
        'Type of Interaction',
        'N.of interactions',
        'Interaction strength'
    ])

    # Iterate through rows of df and store bacteriocin contact residues as a separate index associated
    # with a given receptor contact residue
    for index, row in df.iterrows():
        value = row[0]
        if isinstance(value, str):
            if "pdb" in value:
                new_row = pd.DataFrame([{
                    'PDB file number': value,
                    'Receptor residue and number': '-',
                    'Bacteriocin residue and number': '-',
                    'Type of Interaction': '-',
                    'N.of interactions': '-',
                    'Interaction strength': '-'
                }])
                output_df = pd.concat([output_df, new_row], ignore_index=True)
            elif value in ("Hydrophobic", "PiPi", "CationPi", "Ionic"):
                type_of_interaction = value
            else:
                receptor_residue_to_add = value
        elif isinstance(value, list):  # Check if the entry is a list (bacteriocin contacts)
            # Create a new DataFrame for the new row
            new_row = pd.DataFrame([{
                'PDB file number': '-',
                'Receptor residue and number': receptor_residue_to_add,
                'Bacteriocin residue and number': value[0],
                'Type of Interaction': type_of_interaction,
                'N.of interactions': value[1],
                'Interaction strength': value[2]
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
