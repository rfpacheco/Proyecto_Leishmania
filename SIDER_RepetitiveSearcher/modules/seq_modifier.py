import csv
import subprocess

from modules.files_manager import csv_creator

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 2) Utiliza un archivo CSV, del cual lee las secuencias y las extiende hasta los 1000 nt, creando un archivo CSV resultante con esos datos


def specific_sequence_1000nt(path_input, chromosome_ID, main_folder_path):
    """
    This function will expand the selected sequence to 1000 nt and will create a CSV file with it.

    It uses the function :func:`~modules.file_manager.csv_creator`

    :param path_input: Path where the .csv file we'll use to read and filter data is.
    :type path_input: string

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*.
    :type chromosome_ID: string

    :param main_folder_path: Path where the results will be placed. It will create a subfolder with the chromosome_ID name + "_1000nt.csv".
    :type main_folder_path: string

    :return: a csv file with the sufix "_1000nt.csv"
    :rtype: CSV file
    """

    chrX_1000nt = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")

        for row in reader:

            # -----------------------------------------------------------------------------
            if "plus" in row[14]:
                seq_length = int(row[11]) - int(row[10])  # Way better than using row[3] "Alignment length", because it can stay the same even though the rest DID change with "blastdbcmd".
                # That's why we use row[10] "Start of alignment in subject" and row[11] "End of alignment in subject".
                number_add_length = int((1000 - seq_length) / 2)  # Now we need the coordinates to expand till 1000 nt.
                new_start = int(row[10]) - number_add_length  # New coordinates for Start
                new_end = int(row[11]) + number_add_length  # New coordinates for End

                # Now we select the sequence with new coordinates with "blastdbcmd"
                seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry "
                                              + row[1] + " -range " + str(new_start) + "-" + str(new_end)
                                              + " -strand plus -outfmt %s",
                                              shell=True,
                                              universal_newlines=True)  # Very important to use "subprocess.check_output" so we can get the output
                seq = seq.strip()  # Remove EoL characters

                # We create a custom made CSV row with the new info
                new_row = [row[0], row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                # And we add if to our chrX_1000nt list
                chrX_1000nt.append(new_row)

            # -----------------------------------------------------------------------------
            # Now we do the same with the "minus" strand
            elif "minus" in row[14]:
                seq_length = int(row[10]) - int(row[11])  # This time the rest is inverted compared to "plus". We can as well make abs() instead.

                number_add_length = int((1000 - seq_length)/ 2)
                new_start = int(row[10]) + number_add_length  # Upside down because its "minus"
                new_end = int(row[11]) - number_add_length  # Upside down because its "minus"

                seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " 
                                              + row[1] + " -range " + str(new_end) + "-" + str(new_start) 
                                              + " -strand minus -outfmt %s",
                                              shell=True,
                                              universal_newlines=True)
                seq = seq.strip()

                new_row = [row[0], row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                chrX_1000nt.append(new_row)

                # -----------------------------------------------------------------------------

    # We create a CSV file called 1000nt.csv
    writing_path_input = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_1000nt.csv"
    csv_creator(writing_path_input, chrX_1000nt)

    return (writing_path_input)  # Important to know the path to this file.
