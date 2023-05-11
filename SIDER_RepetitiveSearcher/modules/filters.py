import csv

from Bio import SeqIO

from modules.files_manager import csv_creator
from modules.duplicates import genome_duplicate_filter
from modules.overlap import genome_solap_main

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def chromosome_filter(path_input, name):
    """
    With this filter we obtain the labels/titles for each cromosome of our file, e.g., in **Leishmania** case we'll obtain the labels "LinJ.01", "LinJ.02", etc.
    This filter reads all the sequences in a FASTA file. Then with a prefix ``name``, it adds a numbering in the format ``.XX``, being X the numbers in order for each sequences it finds.


    For example, for:

    - ``name = "LinJ"``
    - A fasta file of 36 sequences:

    The output would be:

    .. code-block:: bash

       ['LinJ.01', 'LinJ.02', 'LinJ.03', 'LinJ.04', 'LinJ.05', 'LinJ.06', 'LinJ.07', 'LinJ.08', 'LinJ.09', 'LinJ.10', 'LinJ.11', 'LinJ.12', 'LinJ.13', 'LinJ.14', 'LinJ.15', 'LinJ.16', 'LinJ.17', 'LinJ.18', 'LinJ.19', 'LinJ.20', 'LinJ.21', 'LinJ.22', 'LinJ.23', 'LinJ.24', 'LinJ.25', 'LinJ.26', 'LinJ.27', 'LinJ.28', 'LinJ.29', 'LinJ.30', 'LinJ.31', 'LinJ.32', 'LinJ.33', 'LinJ.34', 'LinJ.35', 'LinJ.36']

    The objective of this function is to be able to correctly name the files since programming languages do not usually admit a non-string format, numbers that start at 0, in this way we can automate their correct labeling, especially for numbers from 01 to 09 .

    :param path_input: Path to the ``.fasta`` file to read.
    :type path_input: string

    :param name: Name to give the results. In **Leishmania**'s case, it's "LinJ".
    :type name: string

    :return: A python list with the chosen labels
    :rtype: Python list
    """
    max_chr = len(list(SeqIO.parse(path_input, "fasta")))  # It reads the FASTA file and gets the total number of chromosomes.
    chromosome_number = []
    main_list = (list(range(1, max_chr + 1)))  # We index correctly with "+ 1" since Python starts everython in 0
    for number in main_list:
        number = str(number)
        if len(number) == 1:
            chromosome_number.append(name + ".0" + number)  # Para que coincida con los CSV
        else:
            chromosome_number.append(name + "." + number)

    return (chromosome_number)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def filter_by_column(path_input, column, size_filter, writing_path_input):
    """
    This function will filter a CSV data depending on ``length`` (if we want to firlter by sequence length) or ``percent`` (if we want to filter by identity percent).

    :param path_input: Path to the CSV file we want to filter data. It's the output file created by :func:`~modules.blaster.blastn_blaster`.
    :type path_input: string

    :param column: Can be ``length`` (if we want to firlter by sequence length) or ``percent`` (if we want to filter by identity percent)
    :type column: string

    :param size_filter: Number to filter dependint of the ``column`` argument.
    :type size_filter: integer

    :param writing_path_input: Path to the CSV file this function will create and save
    :type writing_path_input: string

    :return: A CSV file with the dalta filtered depending on the ``column`` and ``size_filter`` argumetns.
    :rtype: CSV file
    """

    if column == "length":
        column = 3
    elif column == "percent":
        column = 2

    matrix_filter_by_column = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if column == 3:
                if 1000 >= int(row[column]) >= size_filter:  # 1000 is an important number for the Duplication problem (J.M. Requena Rolania, personal communication)
                    matrix_filter_by_column.append(row)
            elif column == 2:
                if float(row[column]) >= size_filter:  # Needed to go from string to floar
                    matrix_filter_by_column.append(row)

    csv_creator(writing_path_input, matrix_filter_by_column)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def dash_filter(path_input, writing_path_input):
    """
    With this function will filter all the "-" dashes in the sequences.

    :param path_input: Path to a CSV file we want to filter the data from. The data came from the output of :func:`~modelos.filters.filter_by_column`.
    :type path_input: string

    :param writing_path_input: Path to a CSV file where the filtered data will be written.
    :type writing_path_input: string

    :return: A CSV file with all the "-" dashes filtered.
    :rtype: CSV file
    """
    matrix_dash_filter = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            row[15] = row[15].replace("-", "")
            matrix_dash_filter.append(row)

    csv_creator(writing_path_input, matrix_dash_filter)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def global_filters_main(path_input, writing_path_input, genome_fasta, naming_short, max_diff):
    """
    This function mixes every other filter made. Each one writes a CSV which is overwritten every time till the final step.

    :param path_input: Path to the CSV file we want to filter data. It's the output file created by :func:`~modules.blaster.blastn_blaster`.
    :type path_input: string

    :param writing_path_input: Path where the CSV file will be saved.
    :type writing_path_input: string

    :param genome_fasta: Path to our whole genome sequence in FASTA format.
    :type genome_fasta: string

    :param naming_short: Label needed to read the ID of each cromosome in the .csv file. In the case of **L. infantum** for example, would be *LinJ* since the .csv file IDs are *LinJ.XX*.
    :type naming_short: string

    :param max_diff: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type max_diff: intenger

    :return:
    :rtype:
    """

    column = "length"
    size_filter = 100

    # ArithmeticErrorThis will take name "X" CSV file "_BLAST_MAIN.csv" and it will overwrite it with the same name "X"
    filter_by_column(path_input, column, size_filter, writing_path_input)

    path_input = writing_path_input  # This way we tell the program the input file "path_input" is the same as the output file of "filter_by_column". Tbh it's not needed, but its like to improve the understanding.

    # This will take name "X" CSV file "_BLAST_MAIN.csv" and it will overwrite it with the same name "X". So path_input is the same as writing_path_input
    dash_filter(path_input, writing_path_input)

    # Remember "path_input" is the same as "writing_path_input"
    genome_duplicate_filter(genome_fasta, naming_short, path_input, max_diff, writing_path_input)

    genome_solap_main(genome_fasta, naming_short, path_input, max_diff, writing_path_input)
    # En este ultima ya imprime los resultados finales

# global_filters_main(path_input, writing_path_input, genome_fasta, naming_short, max_diff):

    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, Recordar la extension .csv
    # Arg 2: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    # Arg 3: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 4: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE