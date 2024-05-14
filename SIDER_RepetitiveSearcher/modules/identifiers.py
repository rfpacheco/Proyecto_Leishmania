import os
import csv

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

    # -----------------------------------------------------------------------------
    nucleotides1000_directory = specific_sequence_1000nt(data_input, chromosome_ID, main_folder_path, genome_fasta)  # Extend sequence to 1000 nt.

    fasta_creator_output = main_folder_path + chromosome_ID + "/" + chromosome_ID + "_1000nt.fasta"
    fasta_creator(nucleotides1000_directory, fasta_creator_output)

    blastn_dic(fasta_creator_output)

    blaster_output = main_folder_path + chromosome_ID + "/" + chromosome_ID + "_1000nt_Blaster.csv"
    blastn_blaster(fasta_creator_output,  # They are the same
                   fasta_creator_output,  # Because we are launching one against each other --> this way we'll get the correc coordinates in the future.
                   blaster_output,
                   85)

    filter_by_column(blaster_output,
                     "length",
                     100,
                     blaster_output)

    # -----------------------------------------------------------------------------
    corrected_sequences = specific_sequence_corrected(blaster_output, nucleotides1000_directory, main_folder_path, chromosome_ID, genome_fasta)

    # -----------------------------------------------------------------------------
    # This module doesn't work --> need to be redone
    # subfamilies_file_path_writing = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_Subfamily.csv"
    # subfamily_sorter(blaster_output, corrected_sequences, subfamilies_file_path_writing)

    # -----------------------------------------------------------------------------
    # pdb.set_trace()
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