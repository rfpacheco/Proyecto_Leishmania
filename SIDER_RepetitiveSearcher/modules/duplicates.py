import csv

# from modules.filters import chromosome_filter  # Don't call -> circular import
from modules.files_manager import csv_creator


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense, max_diff):
    """
    This one is a filter for duplications in the whole genome. First we check the neighborhood sequences inside one chromosome (by coordinates) and search for equal sequences insice that neighborhood.

    :param genome_fasta: Path to our whole genome sequence in FASTA format.
    :type genome_fasta: string

    :param naming_short: Label needed to read the ID of each cromosome in the .csv file. In the case of **L. infantum** for example, would be *LinJ* since the .csv file IDs are *LinJ.XX*.
    :type naming_short: string

    :param path_input: Path to the CSV file we want to filter data. It's the output file created by :func:`~modules.filters.dash_filter` and given by :func:`~modules.filters.genome_duplicate_filter`.
    :type path_input:

    :param DNA_sense: Be it ``plus`` be it ``minus``, it's used to differentiate the DNA strand.
    :type DNA_sense: string

    :param max_diff: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type max_diff: integer

    :return: All the rows from the CSV file without duplications in a Python matrix.
    :rtype: A Python matrix
    """
    from modules.filters import chromosome_filter  # Delayed import --> to break the ciruclar import. Need to be at the start of function.

    matrix_all_genome = []
    chromosome_number = chromosome_filter(genome_fasta, naming_short)  # It obtains a list with the chromosomes IDs, e.g., ["LinJ.01", "LinJ.02", ...]

    for chromosome in chromosome_number:  # For each chromosome, e.g., it will start only with "LinJ.01"
        location_start = []  # It will save ALL "start" (row[10]) coordinates for a chromosome ID (e.g., "LinJ.01") and a DNA sense (e.g., "plus").
        with open(path_input, "r") as main_file:  # Opens CSV file "BLAST_MAIN.csv" --> after getting the correct coordinates, it has the results from BLASTN of the correct coordinates agains the whole genome.
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:  # Chromosome filter, e.g., it checks if "LinJ.01" is in row[1].
                    if DNA_sense in row[14]:  # Filter correct DNA sense ("plus" or "minus")
                        location_start.append(int(row[10]))  # It saves the "Start of alignment in query" from the CSV file to "location_start" list.

            # -----------------------------------------------------------------------------
            matrix_filter = []
            position_global = []  # VERY IMPORTANT! With this we prevent repeated locations. This one will reset for each chromosome ID in "chromosome_number".
            with open(path_input, "r") as main_file:  # we read the CSV "BLAST_MAIN.csv" again
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    if chromosome in row[1]:  # Chromosome filter
                        if DNA_sense in row[14] and int(row[10]) not in position_global:  # 1) Checks for DNA sense. 2) Checks if "start" is in "position_global" to avoid repeated locations.
                            position_rec = []  # With a specific chromosome ID (e.g., "LinJ.01"), DNA sense, and not in "position_global" --> Saves "positions" from "location_start".

                            for position in location_start:  # "location_start" --> ALL "start" coordinates for a chromosome ID (e.g., "LinJ.01") and a DNA sense (e.g., "plus") are saved here.
                                if abs(int(row[10]) - position) < max_diff:  # Compare row[10] from "BLAST_MAIN.csv" with each "position" in "location_start"
                                    # In this part we make sure we are NEAR our position in the genome (e.g. with max_diff = 600): 
                                        # If position is 35000 and row[10] is 10000, then abs(10000-35000)=25000 > 600, so we are FAR AWAY.
                                        # If the "position" is 35000 and row[10] is 35400, then abs(35400-35000)=400 < 600, so we are near --> It'll try to save the "position".
                                    if position not in position_rec:  # Checks for repeated locations inside "position_rec", which resets for each chromosome.
                                        # We add all non-repeated positions near this row[10]
                                        position_rec.append(position)
                                        position_global.append(position)

                            # For a chromosome we've got "position_rec" and "position_global".
                            DNAseq_filter = []
                            with open(path_input, "r") as main_file:  # Need to open it again to start reading from the first row
                                reader = csv.reader(main_file, delimiter=",")
                                for row in reader:
                                    if chromosome in row[1]:  # Chromosome filter
                                        if int(row[10]) in position_rec:  # if it's in "position_rec", then we check for duplications, since they are more or less near each other.

                                            if row[15] in DNAseq_filter:
                                                continue  # If the seq is inside our "DNAseq_filter", it will skip the rest of the code and jump back to the last loop "for". This way we avoid adding duplications.
                                            else:  # If it's not inside "DNAseq_filter", we add it and to "matrix_filter" as well
                                                DNAseq_filter.append(row[15])
                                                matrix_filter.append(row)
        matrix_all_genome += matrix_filter

    return (matrix_all_genome)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_duplicate_filter(genome_fasta, naming_short, path_input, max_diff, writing_path_input):  # Todo STRING menos max_diff
    """
    This functions just calls :func:`~modules.duplicates.genome_pre_duplicate_filter` for the *plus* and *minus* DNA strand, unites them in a matrix, and creates a CSV file with it.

    :param genome_fasta: Path to our whole genome sequence in FASTA format.
    :type genome_fasta: string

    :param naming_short: Label needed to read the ID of each cromosome in the .csv file. In the case of **L. infantum** for example, would be *LinJ* since the .csv file IDs are *LinJ.XX*.
    :type naming_short: string

    :param path_input: Path to the CSV file we want to filter data. It's the output file created by :func:`~modules.filters.dash_filter`.
    :type path_input:

    :param max_diff: Maximun proxomity value for the different sequences when they have to be grouped. **Important**.
    :type max_diff: integer

    :param writing_path_input: Path where the CSV file will be saved.
    :type writing_path_input: string

    :return: A CSV file with all the data filtered without duplications.
    :rtype: CSV file
    """
    print("\n", "=" * 50, "\nFiltering duplicates proceeding:\n", "=" * 50, sep="")

    DNA_sense = ["plus", "minus"]  # To differentiate between "+" and "-" strand

    # We call first the "plus" strand
    matrix_main1 = genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense[0], max_diff)

    # And then the "minus" strand
    matrix_main2 = genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense[1], max_diff)

    # Then we add one matrix after the other
    matrix_main = matrix_main1 + matrix_main2

    csv_creator(writing_path_input, matrix_main)
