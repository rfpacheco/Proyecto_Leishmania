import os
import time
import pandas as pd

from modules.files_manager import fasta_creator, columns_to_numeric
# from modules.blaster import blastn_dic  # IMPORTANT -> Since "blaster.py" is importing "identifiers.py" I can't make "identifiers.py" import "blaster.py" --> ERROR: CIRCULAR IMPORT
from modules.seq_modifier import specific_sequence_1000nt, specific_sequence_corrected
from modules.filters import global_filters_main
from modules.bedops import bedops_main  # New module 19/04/2024
from modules.stopping import stopping_main  # New module
# from modules.subfamilies_finder import subfamily_sorter  # Needs to be modified


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def genome_specific_chromosome_main(data_input, chromosome_ID, main_folder_path, genome_fasta, identity_1, identity_2, run_phase):
    """
    Main program, which calls the iterative blaster and all the filters needed.

    :param data_input: data frame. Initially is a pandas gropy object wich was output from :func:`~modules.blaster.blastn_blaster` alone to obtain the initial data.
    :type data_input: pandas DataFrame

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*, which were present more then once in the ``path_input``. See :func:`~modules.blaster.repetitive_blaster` for more information.
    :type chromosome_ID: member of a python list

    :param main_folder_path: Path where the results will be placed.
    :type main_folder_path: string

    :param genome_fasta: Path to our whole genome sequence in FASTA format.
    :type genome_fasta: string

    :return: A CSV file with all the data filtered
    :rtype: CSV file
    """
    from modules.blaster import blastn_dic, blastn_blaster  # Delayed import --> to break the ciruclar import. Need to be at the start of function.

    chromosme_folder_path = os.path.join(main_folder_path, chromosome_ID)  # For "chromosome_ID" it creates a folder in the main folder.
    os.makedirs(chromosme_folder_path, exist_ok=True)  # Folder
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    sequences_1000 = specific_sequence_1000nt(data_input, chromosome_ID, main_folder_path, genome_fasta)  # Extend sequence to 1000 nt. Saved into a pandas Data Frame. It modifies the original data_input
    sequences_1000_fasta_path = os.path.join(chromosme_folder_path, chromosome_ID + "_1000nt.fasta")  # Path to the output FASTA file
    toc = time.perf_counter()
    print("")
    print("\t\t2.1. Sequence extension to 1000nt:\n",
          f"\t\t\t- Data row length: {sequences_1000.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    fasta_creator(sequences_1000, sequences_1000_fasta_path)
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.2. Fasta 1000nt file creation:\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    blastn_dic(sequences_1000_fasta_path, sequences_1000_fasta_path)
    tic = time.perf_counter()
    second_blaster = blastn_blaster(sequences_1000_fasta_path, sequences_1000_fasta_path, identity_2)
    second_blaster = columns_to_numeric(second_blaster, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    second_blaster_filtered = second_blaster[second_blaster["length"].astype(int) > 100]  # Filter by length
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.3. Second BLASTn 1000vs1000:\n",
          f"\t\t\t- Data row length: {second_blaster.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Create new folders
    folder_corrected_sequences_path = os.path.join(chromosme_folder_path, "corrected_sequences")
    os.makedirs(folder_corrected_sequences_path, exist_ok=True)
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    corrected_sequences = specific_sequence_corrected(data_input=second_blaster_filtered,
                                                      nucleotides1000_df=sequences_1000,
                                                      first_data_input = data_input,  # this list is modified by `specific_sequence_1000nt`, so most of the data has 1000nt long
                                                      main_folder_path=folder_corrected_sequences_path,
                                                      genome_fasta=genome_fasta,
                                                      chromosome_ID=chromosome_ID,
                                                      run_phase=run_phase)
    corrected_sequences.sort_values(by=["sstrand", "sstart"], inplace=True)  # Sort the data frame inplace
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.4. Corrected sequences:\n",
            f"\t\t\t- Data row length: {corrected_sequences.shape[0]}\n",
            f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Save the corrected sequences to a CSV filecorrected_sequences_path
    corrected_sequences_path = os.path.join(folder_corrected_sequences_path, f"{run_phase}_corrected.csv")
    corrected_sequences.to_csv(corrected_sequences_path, index=False, header=True, sep=",")  # Save the data frame to a CSV file
    # -----------------------------------------------------------------------------
    # Make the checks
    corrected_sequences_df1_path = os.path.join(folder_corrected_sequences_path, f"{run_phase - 1}_corrected.csv")
    corrected_sequences_df2_path = corrected_sequences_path

    if run_phase > 1:  # because run 0 doesn't exist
        corrected_sequences_df1 = pd.read_csv(corrected_sequences_df1_path, sep=",", header=0)
        corrected_sequences_df2 = pd.read_csv(corrected_sequences_df2_path, sep=",", header=0)
        stop_data = stopping_main(data_df1=corrected_sequences_df1, data_df2=corrected_sequences_df2)  # Check if the data is the same as the last one. TRUE or FALSE for this chromosome.
    else:
        stop_data = False
    # -----------------------------------------------------------------------------
    second_fasta_creator_path = os.path.join(chromosme_folder_path, chromosome_ID + "_Corrected.fasta")
    tic = time.perf_counter()
    fasta_creator(corrected_sequences, second_fasta_creator_path)
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.5. Corrected sequences fasta file creation:\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    third_blaster = blastn_blaster(query_path=second_fasta_creator_path,
                                    dict_path=genome_fasta,
                                    perc_identity=identity_1)
    third_blaster = columns_to_numeric(third_blaster, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.6. BLASTn against genome:\n",
          f"\t\t\t- Data row length: {third_blaster.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    filtered_data = global_filters_main(data_input=third_blaster,
                                        genome_fasta=genome_fasta,
                                        writing_path=chromosme_folder_path)
    # filtered_data = columns_to_numeric(filtered_data, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    toc = time.perf_counter()
    print("")
    print("\t\t2.7. Filtering BLASTn against genome:\n",
          f"\t\t\t- Data row length: {filtered_data.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    return filtered_data, stop_data  # Returns the data frame