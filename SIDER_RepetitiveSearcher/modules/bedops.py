# Modules needed
import numpy as np
import pandas as pd
import subprocess
import os

from modules.files_manager import columns_to_numeric

def get_data_sequence(data, strand, genome_fasta):
    """
    This function gets the sequence of the data from the fasta file. It will keep the Chromosome ID, start coordinate, end coordinate and strand.

    :param data: A pandas data frame with the data read of the BED files.
    :type data: pandas.core.frame.DataFrame

    :param strand: The strand of the sequence. It can be "plus" or "minus".
    :type strand: str

    :param genome_fasta: Path to the whole genome sequence in FASTA format.
    :type genome_fasta: string
    """
    sequences = []
    for _, row in data.iterrows():
        sseqid = row["sseqid"]
        start = row["sstart"]
        end = row["send"]
        cmd = [
            "blastdbcmd",
            "-db", genome_fasta,
            "-entry", sseqid,
            "-range", f"{start}-{end}",
            "-strand", strand,
            "-outfmt", "%s"
        ]

        sequence = subprocess.check_output(cmd, universal_newlines=True).replace('\n', '')

        sequences.append({
            "sseqid": sseqid,
            "sstart": start,
            "send": end,
            "sstrand": strand,
            "sseq": sequence
        })

    sequences_df = pd.DataFrame(sequences)
    return sequences_df

def bedops_main(data_input, genome_fasta, writing_path_input):
    """
    This function will implement BEDOPS to filter duplicates and overlaps in a CSV file.

    :param path_input: Path to the CSV file we want to filter data.
    :type path_input: string

    :param genome_fasta: Path to the whole genome sequence in FASTA format.
    :type genome_fasta: string
    """
    # -----------------------------------------------------------------------------
    # 1) Filter and sort data
    # -----------------------------------------------------------------------------
    columns_ids = data_input.columns  # gets the columns names
    column_length = len(columns_ids)  # gets the number of columns
    df_plus = data_input[data_input["sstrand"] == "plus"]  # filters the "+" strand
    df_minus = data_input[data_input["sstrand"] == "minus"]  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=["sseqid", "sstart"])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=["sseqid", "sstart"])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 2) BEDOPS files creation:
    # -----------------------------------------------------------------------------
    #  row[1] == Chromosome ID, row[10] == Start coordinate, row[11] == End coordinate
    plus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_plus.bed")
    minus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_minus.bed")
    
    df_plus[["sseqid", "sstart", "send"]].to_csv(plus_path, sep="\t", header=False, index=False)  # creates a BED file for the "+" strand
    # In the minus strand its important to change the order of the coordinates, because "bedops" reads them like "plus" strand.
    # If not, it will not merge them.
    df_minus[["sseqid", "send", "sstart"]].to_csv(minus_path, sep="\t", header=False, index=False)  # creates a BED file for the "-" strand

    # -----------------------------------------------------------------------------
    # 3) BEDOPS function call with subprocess.
    # -----------------------------------------------------------------------------
    # BEDOPS call to "plus.bed" and "minus.bed" files
    # Using subprocess to call BEDOPS.
    # Using subprocess.check_output() to get the output from the command.
    # shell=True to have the input of .check_output() as a string.
    # universal_newlines=True to have the EoL character as "\n" (Unix-like systems).
    # The output will be a variable of strings.
    df_plus_bedops = subprocess.check_output(f"bedops --merge {plus_path}", shell=True, universal_newlines=True)  # merges the "+" strand BED file
    df_minus_bedops = subprocess.check_output(f"bedops --merge {minus_path}", shell=True, universal_newlines=True)  # merges the "-" strand BED file

     # Now let's transform then into Data Frames
    df_plus_bedops = pd.DataFrame([x.split("\t") for x in df_plus_bedops.split("\n") if x],
                                    columns=["sseqid", "sstart", "send"])  # transforms the "+" strand BEDOPS output into a Data Frame
    df_minus_bedops = pd.DataFrame([x.split("\t") for x in df_minus_bedops.split("\n") if x],
                                    columns=["sseqid", "sstart", "send"])  # transforms the "-" strand BEDOPS output into a Data Frame
    # -----------------------------------------------------------------------------
    # 4) Call `blastdbcmd` to get the sequences with the function get_data_sequence()
    # -----------------------------------------------------------------------------
    if df_plus_bedops.empty:  # In case the original data is empty, the code needs to keep going
        df_plus_bedops_wseq = pd.DataFrame(columns=columns_ids)  # creates an empty Data Frame with 5 columns
    else:  # If the original data is not empty, tit uses get_data_sequence
        df_plus_bedops_wseq = get_data_sequence(df_plus_bedops, "plus", genome_fasta)

    # The same for the minus strand:
    if df_minus_bedops.empty:
        df_minus_bedops_wseq = pd.DataFrame(columns=columns_ids)
    else:   
        df_minus_bedops_wseq = get_data_sequence(df_minus_bedops, "minus", genome_fasta)

 
    # Let's reorderthe `df_minus_bedps_wseq` data frame:
    df_minus_bedops_wseq[["sstart", "send"]] = df_minus_bedops_wseq[["send", "sstart"]].copy()  # swap only values
    # -----------------------------------------------------------------------------
    # 5) Processing data
    # -----------------------------------------------------------------------------
    # Join both data frames
    all_data = pd.concat([df_plus_bedops_wseq, df_minus_bedops_wseq], ignore_index=True)  # joins both Data Frames

    # Adding sequence length to the DataFrame:
    new_column = [len(x) for x in all_data.loc[:,"sseq"]]  # creates a list with the length of each sequence
    all_data.insert(1, "length", new_column)  # inserts the new column with the sequence length. Column index are shifted.


    # -----------------------------------------------------------------------------
    # 6) Correctly modeling the output Data Frame to 15 columns and output as CSV file.
    # -----------------------------------------------------------------------------
    new_data = pd.DataFrame(index=range(all_data.shape[0]), columns=columns_ids)  # creates a new Data Frame with 15 columns. The rows depends on the .shape[0]

    new_data.loc[:,["sseqid", "length", "sstart", "send", "sstrand", "sseq"]] = all_data.loc[:,["sseqid", "length", "sstart", "send", "sstrand", "sseq"]].copy()
    new_data = columns_to_numeric(new_data, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    
    return new_data  # returns the new Data Frame
