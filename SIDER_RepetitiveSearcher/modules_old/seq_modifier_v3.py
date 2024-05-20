import pandas as pd
import subprocess
import re

from modules.files_manager import csv_creator
from modules.bedops import bedops_main, bedops_second

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def specific_sequence_1000nt(data_input, chromosome_ID, main_folder_path, genome_fasta):
    """
    This function will expand the selected sequences to 1000 nt and will create a CSV file with it. For this, it will use ``blastdbcmd`` from BLAST `Command Line Application User Manual`_, which will get the sequences in the **whole fasta genome** file so every sequence returned will be "plus".


    It uses the function :func:`~modules.file_manager.csv_creator`

    :param path_input: Path to the main CSV file where data will be filtered. In this case is the output of expanded sequences to 1000nt from :func:`~modules.identifiers.specific_sequence_extractor`.
    :type path_input: string

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*, which were present more then once in the ``path_input``. See :func:`~modules.blaster.repetitive_blaster` for more information.
    :type chromosome_ID: member of a python list

    :param main_folder_path: Path where the results will be placed. It will create a subfolder with the chromosome_ID name + "_1000nt.csv".
    :type main_folder_path: string

    :return: a CSV file with the sufix "_1000nt.csv"
    :rtype: CSV file
    """
    # -----------------------------------------------------------------------------
    for index, element in data_input.iterrows():
        if "plus" in element["strand"]:
            lower_coor = element["start"].astype(int)  # We get the start of the sequence
            upper_coor = element["send"].astype(int)  # We get the end of the sequence
        else:  # If it's the "-" strand
            lower_coor = element["send"].astype(int)
            upper_coor = element["start"].astype(int)
        
        subject_length = upper_coor - lower_coor + 1
        if subject_length < 1000:  # If the sequence is less than 1000nt, we'll expand it.
            leftover_length = 1000 - subject_length  # We get the difference between 1000 and the length of the sequence
            leftover_length_halved = leftover_length / 2  # We divide it by 2 to get the half of the difference
            lower_coor = lower_coor - leftover_length_halved  # We subtract the half of the difference to the start
            upper_coor = upper_coor + leftover_length_halved  # We add the half of the difference to the end

            # Check if the coordinates are not out of the genome
            if lower_coor < 0:
                lower_coor = 1
                upper_coor += leftover_length_halved  # Since `new_start` is 1, we add the half of the difference to the end
            
            # We get the sequence from the whole genome
            cmd = f"blastdbcmd -db {genome_fasta} -entry {element['seqid']} -range {lower_coor}-{upper_coor} -strand {element["strand"]} -outfmt %s"
            seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

            # Now with the data let's modify the data frame
            data_input.loc[index, "pident"] = ""
            data_input.loc[index, "length"] = ""
            data_input.loc[index, "qsstart"] = ""
            data_input.loc[index, "qend"] = ""

            if "plus" in element["strand"]:
                data_input.loc[index, "sstart"] = lower_coor
                data_input.loc[index, "send"] = upper_coor
            else:
                data_input.loc[index, "sstart"] = upper_coor
                data_input.loc[index, "send"] = lower_coor

            data_input.loc[index, "evalue"] = ""
            data_input.loc[index, "bitscore"] = ""
            data_input.loc[index, "qlen"] = ""
            data_input.loc[index, "slen"] = len(seq)
            data_input.loc[index, "sseq"] = seq
        else:  # If the sequence is already 1000nt we just pass
            pass
        
    return data_input
            
 

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def specific_sequence_corrected(data_input, nucleotides1000_directory, main_folder_path, chromosome_ID, genome_fasta):
    """
    The main use of this function is to get the real "coordinates" of the sequence.


    1. First we read the 1000nt CSV file. And we extract from row[0] (i.e., Query row) every *different* query (not repeated.
    2. Then we open the 1000nt_BLASTER file from the 1000nt against each other. To understand it:

       - In the 1000nt file there are "x" number of rows. Each row will be given an index like "Seq_z_LinJ.01_plus", were "z" is the index.
       - So the first row (row[0]) in the 1000nt file will be given "Seq_1.." as a name. And every row with the "Seq_1.." subject in the 1000nt_BLASTER file will the result from the first row in the 1000nt file (row[0]).
       - There will be a maximum of "x" "Seq_x.." in the 1000nt_BLASTER file.

    3. In the 1000nt_BLASTER we search every result from a specific "Seq_z..". And we get the maximum alignment and it's coordinates.
    4. We'll get the "z" number and search it's corresponding sequence in the 1000nt file. And with that and the coordinates, we use ``blastcmd`` to get the correct coordinates from the 1000nt sequence. We need to be careful if the sequence is "plus" or "minus".

    :param path_input: Path to the CSV file we'll use to filter data. It's the result from a BLASTn made between the expanded 1000nt sequences.
    :type path_input: string

    :param nucleotides1000_directory: *return* result from the function :func:`~specific_sequence_1000nt`
    :type nucleotides1000_directory: string

    :param main_folder_path: Path where we'll save the CSV data file to.
    :type main_folder_path: string

    :param chromosome_ID: Identification of the cromosome, e.g., "LinJ.07"
    :type chromosome_ID: string

    :return: CSV File with the corrected coordinates of our sequences.
    :rtype: CSV file
    """

    # First from the BLASTn (one againts each other), we get the IDs of the sequences.

    data_input_selected = data_input[data_input["qseqid"] != data_input["sseqid"]]  # We filter the data to get only the sequences that are different from each other.
    data_input_rejected = data_input[data_input["qseqid"] == data_input["sseqid"]]  # We filter the data to get only the sequences that are the same.

    selected_grouped = data_input_selected.groupby("sseqid")  # We group the data by the subject ID
    new_df = pd.DataFrame()  # We create an empty data frame to store the data
    for name, group in selected_grouped:
        filtered_df = bedops_second(data_input=group,
                                  writing_path_input=main_folder_path)
        new_df = pd.concat([new_df, filtered_df], ignore_index=True)
    
