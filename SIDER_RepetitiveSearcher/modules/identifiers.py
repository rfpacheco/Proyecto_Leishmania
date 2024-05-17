import os
import time

from modules.files_manager import folder_creator, csv_creator, fasta_creator, csv_mixer
# from modules.blaster import blastn_dic  # IMPORTANT -> Since "blaster.py" is importing "identifiers.py" I can't make "identifiers.py" import "blaster.py" --> ERROR: CIRCULAR IMPORT
from modules.seq_modifier import specific_sequence_1000nt, specific_sequence_corrected
from modules.filters import filter_by_column, global_filters_main
from modules.bedops import bedops_main  # New module 19/04/2024
# from modules.subfamilies_finder import subfamily_sorter  # Needs to be modified


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def genome_specific_chromosome_main(data_input, chromosome_ID, main_folder_path, genome_fasta):
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
    data_input_backup = data_input.copy()  # Backup of the original data frame

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    sequences_1000 = specific_sequence_1000nt(data_input, chromosome_ID, main_folder_path, genome_fasta)  # Extend sequence to 1000 nt. Saved into a pandas Data Frame. It modifies the original data_input
    sequences_1000_fasta_path = os.path.join(chromosme_folder_path, chromosome_ID + "_1000nt.fasta")  # Path to the output FASTA file
    toc = time.perf_counter()
    print(f"==>1000nt data frame row length: {sequences_1000.shape[0]}")
    print(f"==>1000nt data frame sequence creation took {toc - tic:0.2f} seconds")
    
    tic = time.perf_counter()
    fasta_creator(sequences_1000, sequences_1000_fasta_path)
    toc = time.perf_counter()
    print(f"==>Fasta 1000nt file creation took {toc - tic:0.2f} seconds")

    blastn_dic(sequences_1000_fasta_path, sequences_1000_fasta_path)
    tic = time.perf_counter()
    second_blaster = blastn_blaster(sequences_1000_fasta_path, sequences_1000_fasta_path, 100)
    second_blaster_filtered = second_blaster[second_blaster["length"].astype(int) > 100]  # Filter by length
    toc = time.perf_counter()
    print(f"==>Second BLASTn row length: {second_blaster.shape[0]}")
    print(f"==>Second BLASTn step took {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    corrected_sequences = specific_sequence_corrected(data_input=second_blaster_filtered,
                                                      nucleotides1000_df=sequences_1000,
                                                      first_data_input = data_input,  # this list is modified by `specific_sequence_1000nt`, so most of the data has 1000nt long
                                                      main_folder_path=main_folder_path,
                                                      genome_fasta=genome_fasta,
                                                      chromosome_ID=chromosome_ID)
    # corrected_sequences = specific_sequence_corrected(second_blaster_filtered, sequences_1000, chromosme_folder_path, chromosome_ID, genome_fasta)

    # -----------------------------------------------------------------------------
    second_fasta_creator_output = main_folder_path + chromosome_ID + "/" + chromosome_ID + "_Corrected.fasta"
    fasta_creator(corrected_sequences, second_fasta_creator_output)

    second_blaster_output = main_folder_path + chromosome_ID + "/" + chromosome_ID + "_BLAST_MAIN.csv"
    blastn_blaster(second_fasta_creator_output,
                   genome_fasta,
                   second_blaster_output,
                   60)

    # -----------------------------------------------------------------------------
    # # SHOULD NOT BE HERE. IS INSIDE global_filters_main(). NEED TO BE REMOVED
    # bedops_main(second_blaster_output,  # input CSV file
    #             genome_fasta,  # Path to the whole genome sequence in FASTA format
    #             second_blaster_output)  # Output to CSV file


    global_filters_main(second_blaster_output,
                        second_blaster_output,
                        genome_fasta,
                        naming_short,
                        max_diff)

    # -----------------------------------------------------------------------------
    csv_mixer_output = main_folder_path + "MIXER.csv"

    # pdb.set_trace()
    if os.path.isfile(csv_mixer_output) is False:  # When it doesn't exist, we create it
        csv_mixer(path_input, second_blaster_output, csv_mixer_output)  # To mix it
    else:  # If the file already exist (already been created), the path changes to "csv_mixer_outpu"
        csv_mixer(csv_mixer_output, second_blaster_output, csv_mixer_output)