import csv
import subprocess

from modules.files_manager import csv_creator

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


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

                number_add_length = int((1000 - seq_length) / 2)
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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# 3) Corrector se secuencias, obtendra las secuencias originales.


def specific_sequence_corrected(path_input, nucleotides1000_directory, main_folder_path, chromosome_ID):
    """
    .. danger::
       NEED TO MODIFY IT BECAUSE OF SOME ERRORS

    The main use of this function is to get the real "coordinates" of the sequence.

    Remember we've got first, after the initial BLASNn, sequence expanded to 1000nt. The next part of the program was launching them one against each other to get a BLASTn which can provide us useful information about the real coordinates:
    
    Now we've got the 1000 nt sequences with an aproximade location of where the SIDER2 are. And then, we've got the BLASTn (of them against each other) which can provide us with the real coordinates.
    
    With this, this function will get the real coordinates.

    :param path_input:
    :type path_input:

    :param nucleotides1000_directory:
    :type nucleotides1000_directory:

    :param main_folder_path:
    :type main_folder_path:

    :param chromosome_ID:
    :type chromosome_ID:

    :return:
    :rtype:
    """

    # First from the BLASTn (one againts each other), we get the IDs of the sequences.
    names = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if row[0] not in names:  # We select the chromosome ID form row[0], so we can have a no-repeated list of those.
                names.append(row[0])  # And example would be "Seq_2_LinJ.01_plus"

    # -----------------------------------------------------------------------------
    # Now, we'll get from the BLASTn (one againts each other), the best alignment.
    chr_x_corrected = []
    for query in names:  # For each chromosome ID row[0] in "names"
        start = []
        end = []
        diference_end_minu_start = 0  # Check it out later. It's to compare it with "difference"

        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if query in row[1]:  # This is row[1], e.g., "Seq_2_LinJ.01_plus"
                    if row[0] != query:  # This is to remove the sequence that overlaps with itself. So if "LinJ.01" overlaps with "LinJ.01", we don't analyze it further.
                        # We don¡t need to difference between "+" and "-" strands, because after BLASTn all the results behave like "+".
                        difference = int(row[11]) - int(row[10])  # Due to how the code is made. Now row[11] will always be > row[10]
                        # Here we'll iterate until we get the larger "difference", i.e., larger "alignment length".
                        if difference > diference_end_minu_start:  # If it's greater than 0
                            diference_end_minu_start = difference  # We save that difference, i.e., alignment length
                            start = []  # Reset "start"
                            end = []  # Reset "end"
                            start.append(int(row[10]))  # Save "start of alignment in subject"
                            end.append(int(row[11]))  # Save "end of alignment in subject"

        # -----------------------------------------------------------------------------
        # Remember we don't need to difference between "+" and "-".
        # In this part, for a specific "query" we'll have the bigger "alignment length" wit its "start" and "end".

        # With the variable "difference" made, this can be removed
        if len(start) > 0 and len(end) > 0:  # I mean this
            min_start = min(start)  # And this
            max_end = max(end)  # And this last one.

            correct_seq = query  # Changed the name to understand it better for the next part.
            number_for_location = int(correct_seq[4]) - 1  # This way we get the 5th position from "Seq_2_LinJ.01_plus". The 4th position is the "number". For example in "Seq_2..." we'll get the "2". Then we rest 1, because in Python everything starts in 0 and not 1. --> ESTO ESTA MAL, Y SI ES UN  "14" DE DOS CIFRAS

            #  Now we filter "correct_seq" to obtain a number to filter a CSV list of 4 x 1000 without doing BLAST
            rows_by_number = []  # I need this part to know if it's the correct row while doing comparisons. This way I can compare it with "rows_by_number[0]" or "[4]" or "[3]" without going in order. We do this in the CSV 4 x 1000nt before the blaster to themselves.
            with open(nucleotides1000_directory, "r") as main_file:
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    rows_by_number.append(row)  # Here we get all the rows from the CSV

            with open(nucleotides1000_directory, "r") as main_file:
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    if row == rows_by_number[number_for_location]:  # Asi me aseguro que estoy en la adecuada. Quizas es rizar el rizo pero no se me ocurre en el momento un paso mejor --> ESTO ESTA MAL
                        if "plus" in row[14]:

                            x = 1000 - max_end
                            new_start = int(row[10]) + min_start - 1
                            new_end = int(row[11]) - x

                            seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry "
                                                          + row[1] + " -range " + str(new_start) + "-" + str(new_end)
                                                          + " -strand plus -outfmt %s",
                                                          shell=True,
                                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
                            seq = seq.strip()  # Eliminar EoL caracteres

                            new_row = [query, row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                            chr_x_corrected.append(new_row)

                        elif "minus" in row[14]:  # Recordar que por como son las coordenadas de minus, que se definen segun las de la posicion plus en 3' --> 5'
                            x = 1000 - max_end
                            new_start = int(row[10]) - min_start + 1
                            new_end = int(row[11]) + x

                            seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " 
                                                          + row[1] + " -range " + str(new_end) + "-" + str(new_start) 
                                                          + " -strand minus -outfmt %s", 
                                                          shell=True, 
                                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
                            seq = seq.strip()  # Eliminar EoL caracteres

                            new_row = [query, row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                            chr_x_corrected.append(new_row)

        if len(start) == 0 and len(end) == 0:  # para casos en los que solo tenga homologia con el mismo, la secuencia se descarta, pero seria mejor cambiar este codigo para insertarla en los siguientes documentos pero no en la forma de 1000nt
            print("\nALERT: individual " + query + " has no homology with no other seq, so it will not be added to the corrected seqs")

    writing_path_input = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_Corrected.csv"
    csv_creator(writing_path_input, chr_x_corrected)

    return (writing_path_input)

# specific_sequence_corrected(path_input, nucleotides1000_directory, main_folder_path, chromosome_ID)

    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: Resultado del RETURN de la funcion Specific_sequence_1000nt
    # Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa
    # Arg 3: STRING. Identificacion del cromosoma, e.g., "LinJ.07"