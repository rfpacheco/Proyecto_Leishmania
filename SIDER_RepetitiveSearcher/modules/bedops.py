# Modules needed
import numpy as np
import pandas as pd
import subprocess
import os

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
    list = []  # It'll save the rows here.
    for index, row in data.iterrows():  # iteration through data frame
        if strand == "plus":  # If strand is "plus", then `blastdbcmd` will take start as row[1] and end as row[2]
            start = row[1]
            end = row[2]
        else: # If the strand is "minus", then `blastdbcmd` will take start as row[2] and end as row[1]
            start = row[2] 
            end = row[1] 
        # `blastdbcmd` call with subprocess.check_output().
        sequence = subprocess.check_output("blastdbcmd -db " + genome_fasta + 
                                           " -entry " + row[0] + 
                                           " -range " + str(start) + "-" + str(end) + 
                                           " -strand " + strand + 
                                           " -outfmt %s", 
                                           shell=True, universal_newlines=True)
        # Appending subprocess the data to the list
        list.append(row[0] + "," + 
                    str(row[1]) + "," + 
                    str(row[2]) + "," + 
                    strand + "," + 
                    sequence)

    # list values are separated by commas, so we split them and create a data frame
    list_split = [row.split(",") for row in list]  # Splitting the list by commas
    list_split_df = pd.DataFrame(list_split)  # Creating a data frame from the list
    list_split_df[4] = list_split_df[4].str.replace('\n', '')  # Important. It removes the new line character from the sequence.

    return list_split_df  # Returns the data frame

def bedops_main(data_input, genome_fasta, writing_path_input, mode=1):
    """
    This function will implement BEDOPS to filter duplicates and overlaps in a CSV file.

    :param path_input: Path to the CSV file we want to filter data.
    :type path_input: string

    :param genome_fasta: Path to the whole genome sequence in FASTA format.
    :type genome_fasta: string
    """
    # -----------------------------------------------------------------------------
    # 1) Import the data from the CSV file
    # -----------------------------------------------------------------------------
    # df = pd.read_csv(path_input, sep=",", header=None)  # reads the data from the CSV file

    # -----------------------------------------------------------------------------
    # 2) Filter and sort data
    # -----------------------------------------------------------------------------
    columns_ids = data_input.columns  # gets the columns names
    column_length = len(columns_ids)  # gets the number of columns
    df_plus = data_input[data_input["sstrand"] == "plus"]  # filters the "+" strand
    df_minus = data_input[data_input["sstrand"] == "minus"]  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=["sseqid", "sstart"])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=["sseqid", "sstart"])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 3) BEDOPS files creation:
    # -----------------------------------------------------------------------------
    #  row[1] == Chromosome ID, row[10] == Start coordinate, row[11] == End coordinate
    plus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_plus.bed")
    minus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_minus.bed")
    
    df_plus[["sseqid", "sstart", "send"]].to_csv(plus_path, sep="\t", header=False, index=False)  # creates a BED file for the "+" strand
    # In the minus strand its important to change the order of the coordinates, because "bedops" reads them like "plus" strand.
    # If not, it will not merge them.
    df_minus[["sseqid", "send", "sstart"]].to_csv(minus_path, sep="\t", header=False, index=False)  # creates a BED file for the "-" strand

    # -----------------------------------------------------------------------------
    # 4) BEDOPS function call with subprocess.
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
    df_plus_bedops = pd.DataFrame([x.split("\t") for x in df_plus_bedops.split("\n") if x])  # transforms the "+" strand BEDOPS output into a Data Frame
    df_minus_bedops = pd.DataFrame([x.split("\t") for x in df_minus_bedops.split("\n") if x])  # transforms the "-" strand BEDOPS output into a Data Frame  

    # Not, let's reorder the minus strand Data Frame
    if not df_minus_bedops.empty:  # If the "-" strand Data Frame is not empty
        df_minus_bedops = df_minus_bedops[[0, 2, 1]]  # reorders the "-" strand Data Frame
        df_minus_bedops.columns = range(df_minus_bedops.columns.size)  # repairs the column index

    # -----------------------------------------------------------------------------
    # 5) Call `blastdbcmd` to get the sequences with the function get_data_sequence()
    # -----------------------------------------------------------------------------
    if mode == 1:  # If the mode is "full", it will get the sequences
        if df_plus_bedops.empty:  # In case the original data is empty, the code needs to keep going
            df_plus_bedops_wseq = pd.DataFrame(columns=range(5))  # creates an empty Data Frame with 5 columns
        else:  # If the original data is not empty, tit uses get_data_sequence
            df_plus_bedops_wseq = get_data_sequence(df_plus_bedops, "plus", genome_fasta)

        # The same for the minus strand:
        if df_minus_bedops.empty:
            df_minus_bedops_wseq = pd.DataFrame(columns=range(5))
        else:   
            df_minus_bedops_wseq = get_data_sequence(df_minus_bedops, "minus", genome_fasta)

        # -----------------------------------------------------------------------------
        # 6) Processing data
        # -----------------------------------------------------------------------------
        # Join both data frames
        all_data = pd.concat([df_plus_bedops_wseq, df_minus_bedops_wseq], ignore_index=True)  # joins both Data Frames

        # Adding sequence length to the DataFrame:
        new_column = [len(x) for x in all_data[4]]  # creates a list with the length of each sequence
        all_data.insert(1, "New", new_column)  # inserts the new column with the sequence length. Column index are shifted.

        # Repair column index
        all_data.columns =range(all_data.columns.size)  # repairs the column index

        # -----------------------------------------------------------------------------
        # 7) Correctly modeling the output Data Frame to 15 columns and output as CSV file.
        # -----------------------------------------------------------------------------
        new_data = pd.DataFrame(index=range(all_data.shape[0]), columns=range(column_length))  # creates a new Data Frame with 15 columns. The rows depends on the .shape[0]
        new_data.iloc[:, [1, 3, 10, 11, 14, 15]] = all_data.iloc[:, [0, 1, 2, 3, 4, 5]]
        
        return new_data  # returns the new Data Frame
    
    elif mode == 2:  # If the mode is not "full", it will return the BEDOPS output
        # add a new column with the strand values:
        if not df_plus_bedops.empty:
            df_plus_bedops.insert(3, 3, "plus")
        else:
            pass

        if not df_minus_bedops.empty:
            df_minus_bedops.insert(3, 3, "minus")
        else:
            pass
        return pd.concat([df_plus_bedops, df_minus_bedops], ignore_index=True)  # returns the BEDOPS output

def bedops_second(data_input, writing_path_input):
    # -----------------------------------------------------------------------------
    # 2) Filter and sort data
    # -----------------------------------------------------------------------------
    columns_ids = data_input.columns  # gets the columns names
    column_length = len(columns_ids)  # gets the number of columns
    df_plus = data_input[data_input["sstrand"] == "plus"]  # filters the "+" strand
    df_minus = data_input[data_input["sstrand"] == "minus"]  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=["sseqid", "sstart"])  # sorts the "+" strand by the start coordinate

    # Let's correct the "-" strand coordinates  
    minus_start = df_minus["sstart"].copy()
    minus_end = df_minus["send"].copy()
    df_minus.loc[:,"sstart"] = list(minus_end)
    df_minus.loc[:,"send"] = list(minus_start)
    df_minus = df_minus.sort_values(by=["sseqid", "sstart"])  # sorts the "-" strand by the start coordinate
    
    all_data = pd.concat([df_plus, df_minus], ignore_index=True)  # joins both Data Frames

    all_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_all.bed")
    all_data[["sseqid", "sstart", "send"]].to_csv(all_path, sep="\t", header=False, index=False)  # creates a BED file for the "-" strand

    all_bedops = subprocess.check_output(f"bedops --merge {all_path}", shell=True, universal_newlines=True)
    all_bedops = pd.DataFrame([x.split("\t") for x in all_bedops.split("\n") if x])

    return all_bedops
