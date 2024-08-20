"""
File: io_utils.py

Module with a get_filehandle function that can be imported by other modules
"""

import os.path
import sys
import argparse
import re


def get_filehandle(file=None, mode=None):
    """
    filehandle : get_filehandle(infile, "r")
    Takes : 2 arguments file name and mode i.e. what is needed to be done with
    this file. This function opens the file based on the mode passed in
    the argument and returns filehandle.
    @param file: The file to open for the mode
    @param mode: They way to open the file, e.g. reading, writing, etc
    @return: filehandle
    """

    try:
        fobj = open(file, mode, encoding='windows-1252')
        return fobj
    except OSError:
        print(f"Could not open the file: {file} for type '{mode}'",
              file=sys.stderr)
        raise
    except ValueError:
        print(f"Could not open the file: {file} for type '{mode}'",
              file=sys.stderr)
        raise


def extract_pdb_information_from_table(large_table):
    # variable stores the number of pdb files in the contacts csv log
    number_of_files = len(large_table[large_table['PDB file number'].str.contains('pdb', na=False)])
    # create an empty dictionary of lists
    dict_of_list = {val: [] for val in range(number_of_files)}
    # stores the df of all rows with pdb labels in a row
    pdb_fields = large_table[large_table['PDB file number'].str.contains('pdb', na=False)]
    # stores a list containing the row numbers of all rows with pdb labels and the last row in the original df
    search_range = list(pdb_fields.index) + [large_table.index[-1]]

    return dict_of_list, search_range


def find_sequence_number(residue):
    query = re.search(r'\d+$', residue)
    return int(query.group()) if query else float('inf')


def get_cli_args():
    parser = argparse.ArgumentParser(description='Provide path+name of YASARA log file data')

    parser.add_argument('-i', '--input', type=str, help='provide the YASARA log file name', required=True)
    return parser.parse_args()
