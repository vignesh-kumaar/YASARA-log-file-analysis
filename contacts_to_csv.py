"""
File: contacts_to_csv.py

This program takes YASARA log file data (contacts.docx) and outputs necessary
information into a csv file

example command: python3 contacts_to_csv.py -i input_files/
python3 contacts_to_csv.py -i ecoli_bamA_mosselii_llpA1/analysis/
"""

from docx import Document
import pandas as pd
from modules import io_utils


def main():
    """Business Logic"""

    # declaring and initializing data
    args = io_utils.get_cli_args()
    infile = args.input

    # reading log data into python
    df = read_contacts_to_df(infile + 'contacts.docx')

    # writing to output file
    df.to_csv(infile + 'contacts.csv', sep='\t', encoding='utf-8', index=False)


def read_contacts_to_df(file=None):
    """
    Reads required information into an empty list and converts it to a dataframe
    :param file: path to the contacts file + filename
    :return: dataframe with required information from contacts.docx
    """

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
            elif words[0] == "Contacts":
                data_array.append([words[3] + words[5]] + words[7:8])

    # Convert list into df for convenience in the next step
    df = pd.DataFrame(data_array)
    df = _format_df(df)
    return df


def _format_df(df):
    """
    Formats the dataframe to distribute the information over distinct columns
    :param df: dataframe with information from contacts.docx
    :return: dataframe
    """

    # Step 2
    # Initialize variable and new df for the next step
    receptor_residue_to_add = None
    output_df = pd.DataFrame(columns=[
        'PDB file number',
        'Receptor residue and number',
        'Bacteriocin residue and number',
        'Strength of contacts'
    ])

    # Iterate through rows of df and store bacteriocin contact residues as a separate index
    # associated with a given receptor contact residue
    for _, row in df.iterrows():
        value = row[0]
        if isinstance(value, str):
            if "pdb" in value:
                new_row = pd.DataFrame([{
                    'PDB file number': value,
                    'Receptor residue and number': '-',
                    'Bacteriocin residue and number': '-',
                    'Strength of contacts': '-'
                }])
                output_df = pd.concat([output_df, new_row], ignore_index=True)
            else:
                receptor_residue_to_add = value
        elif isinstance(value, list):  # Check if the entry is a list (bacteriocin contacts)
            # Create a new DataFrame for the new row
            new_row = pd.DataFrame([{
                'PDB file number': '-',
                'Receptor residue and number': receptor_residue_to_add,
                'Bacteriocin residue and number': value[0],
                'Strength of contacts': value[1]
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
