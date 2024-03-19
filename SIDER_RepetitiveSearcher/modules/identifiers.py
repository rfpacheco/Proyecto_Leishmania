import pdb
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


def specific_sequence_extractor(path_input, chromosome_ID, main_folder_path):
    """
    First it reads a CSV file and takes the rows which include in the second column (the one with the chromosome IDs), the defined chromosome.
    Then, it creates a .csv file  with all the selectes rows in a subfolder in a specified directory using :func:`~modules.files_manager.folder_creator` and :func:`~modules.files_manager.csv_creator`.

    .. admonition:: Example of use

       If ``chromosome_ID = LinJ.02``, it will take all rows in the **.csv** where the second column contains ``LinJ.02``

    :param path_input: Path to the main CSV file where data will be filtered. Initially is a CSV file wich was output from :func:`~modules.blaster.blastn_blaster` alone to obtain the initial data.
    :type path_input: string

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*, which were present more then once in the ``path_input``. See :func:`~modules.blaster.repetitive_blaster` for more information.
    :type chromosome_ID: member of a python list

    :param main_folder_path: Path where the results will be placed. It will create a subfolder with the ``chromosome_ID`` name.
    :type main_folder_path: string

    :return: A CSV file with the selected "chromosome_ID" rows.
    :rtype: CSV file
    """

    # From the file we get all the rows with a specific "chromosome_ID" in row[1]
    chr_x_seqs = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if chromosome_ID in row[1]:
                chr_x_seqs.append(row)

    # Writing output for a folder creation
    folder_path = main_folder_path + chromosome_ID
    folder_creator(folder_path)

    # In the folder, writing ouput for the CSV file
    writing_path_input = main_folder_path + chromosome_ID + "/" + chromosome_ID + ".csv"
    csv_creator(writing_path_input, chr_x_seqs)

    # We return the folder path and the CSV path
    return (folder_path, writing_path_input)  # I think "folder path" is useless here in the return

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_specific_chromosome_main(path_input, chromosome_ID, main_folder_path, genome_fasta, naming_short, max_diff):
    """
    Main program, which calls the iterative blaster and all the filters needed.

    :param path_input: Path to the main CSV file where data will be filtered. Initially is a CSV file wich was output from :func:`~modules.blaster.blastn_blaster` alone to obtain the initial data.
    :type path_input: string

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*, which were present more then once in the ``path_input``. See :func:`~modules.blaster.repetitive_blaster` for more information.
    :type chromosome_ID: member of a python list

    :param main_folder_path: Path where the results will be placed.
    :type main_folder_path: string

    :param genome_fasta: Path to our whole genome sequence in FASTA format.
    :type genome_fasta: string

    :param naming_short: Label needed to read the ID of each cromosome in the .csv file. In the case of **L. infantum** for example, would be *LinJ* since the .csv file IDs are *LinJ.XX*.
    :type naming_short: string

    :param max_diff: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type max_diff: integer

    :return: A CSV file with all the data filtered
    :rtype: CSV file
    """
    from modules.blaster import blastn_dic, blastn_blaster  # Delayed import --> to break the ciruclar import. Need to be at the start of function.

    new_directories = specific_sequence_extractor(path_input, chromosome_ID, main_folder_path)  # We get a .csv with the specified chromosome_ID.
    # folder_path = new_directories[0]  # Chromosome's directory, i.e., folder_path from return (folder_path, writing_path_input) --> I think this was useless          
    last_output = new_directories[1]  # Chromosome's .csv file inside Chromosome's directory, i.e., writing_path_input from return (folder_path, writing_path_input).

    # -----------------------------------------------------------------------------
    nucleotides1000_directory = specific_sequence_1000nt(last_output, chromosome_ID, main_folder_path, genome_fasta)  # Extend sequence to 1000 nt.

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
    # SHOULD NOT BE HERE. IS INSIDE global_filters_main(). NEED TO BE REMOVED
    bedops_main(second_blaster_output,  # input CSV file
                genome_fasta,  # Path to the whole genome sequence in FASTA format
                second_blaster_output)  # Output to CSV file


    # global_filters_main(second_blaster_output,
    #                     second_blaster_output,
    #                     genome_fasta,
    #                     naming_short,
    #                     max_diff)

    # -----------------------------------------------------------------------------
    csv_mixer_output = main_folder_path + "MIXER.csv"

    # pdb.set_trace()
    if os.path.isfile(csv_mixer_output) is False:  # When it doesn't exist, we create it
        csv_mixer(path_input, second_blaster_output, csv_mixer_output)  # To mix it
    else:  # If the file already exist (already been created), the path changes to "csv_mixer_outpu"
        csv_mixer(csv_mixer_output, second_blaster_output, csv_mixer_output)