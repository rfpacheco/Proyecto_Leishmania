import os
import csv
import pdb  # For debugging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def folder_creator(options):
    """
    Creates a Folder in the current path

    :paran options: Name of the folder
    :type options: string

    :returns: Folder in the current directory
    """
    pdb.set_trace()  # Debugging mark
    options = str(options.file_name)
    if not os.path.exists(options):
        os.mkdir(options)
        print("\nDirectory", options, "created at:\n",
              os.path.abspath(options))
    else:
        print("\nDirectory", options, "already exists")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def csv_creator(writing_path_input, writing_input):
    """
    This function will create .csv files

    :param writing_path_input: Path where the .csv file will be created
    :type writing_path_input: string

    :param writing_input: What info we want to add to the .csv file, normally a Matrix/Array 3D
    :type writing_input: TO BE DONE
    """
    with open(writing_path_input, "w") as OutCSV:
        writer = csv.writer(OutCSV)
        writer.writerows(writing_input)
        print("\nCSV:", writing_path_input, "has been created.")


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def csv_mixer(path_input1, path_input2, writing_path_input):
    """
    Mixes two .csv files. Then, uses :func:`~csv_creator` to generate a .csv with the mixed files.

    :param path_input1: Path to the first .csv to mix. It will take the first rows in the final mixed file.
    :type path_input1: string

    :param path_input2: Path to the second .csv to mix. It will take rows adther the first .csv.
    :type path_input2: string

    :param writing_path_input: path to both .csv files
    :type writing_path_input: string
    """

    csv_mixer_matrix = []
    with open(path_input1, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            csv_mixer_matrix.append(row)
    with open(path_input2, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            csv_mixer_matrix.append(row)

    csv_creator(writing_path_input, csv_mixer_matrix)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def fasta_creator(path_input, fasta_output_path):
    """
    This function will create a FASTA file from the input CSV file. For this it will use the **Biopyton** module.

    :param path_input: Path of the CSV file we want to read to transform it to a FASTA file.
    :type path_input: string

    :param fasta_output_path: Path to where we want to save the FASTA file.
    :type fasta_output_path: string

    :return: All the data from the CSV in a FASTA format.
    :rtype: FASTA File

    .. warning::
       MUST MODIFY IT, I NEED THE NUMBERING TO BE 01, 02, 03, 04...09, 10, 11. And not how I have it now, which is 1, 2, 3, 4, 5...10, 11, 12.
    """
    matrix_fasta_creator = []
    numbering = 0
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            numbering += 1
            rec = SeqRecord(
                Seq(row[15]),
                id="Seq_" + str(numbering) + "_" + row[1] + "_" + row[14],  # Que tenga aqui el sentido es esencial para luego filtrarlos
                description="Leishmania infantum " + row[14]
            )
            matrix_fasta_creator.append(rec)

    SeqIO.write(matrix_fasta_creator, fasta_output_path, "fasta")
    print("\nFasta created at:", fasta_output_path)

# fasta_creator(path_input, fasta_output_path)

    # Arg 0: STRING. Directorio del archivo CSV a leer de donde queremos extraer las secuencias FASTA
    # Arg 1: STRING. Directorio del archivo FASTA que contiene las secuencias del archivo CSV. Recordar terminar en la extension .fasta