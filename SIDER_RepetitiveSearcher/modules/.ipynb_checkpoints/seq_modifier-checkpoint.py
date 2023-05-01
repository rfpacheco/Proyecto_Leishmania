import csv
import subprocess

from modules.files_manager import csv_creator

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 2) Utiliza un archivo CSV, del cual lee las secuencias y las extiende hasta los 1000 nt, creando un archivo CSV resultante con esos datos


def specific_sequence_1000nt(path_input, chromosome_ID, main_folder_path):

    chrX_1000nt = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")

        for row in reader:

            if "plus" in row[14]:
                seq_length = int(row[11]) - int(row[10])  # Way better than using row[3] "Alignment length", because it can stay the same even though the rest DID change with "blastdbcmd".
                # That's why we use row[10] "Start of alignment in subject" and row[11] "End of alignment in subject".
                number_add_length = int((1000 - seq_length) / 2)  # Now we need the coordinates to expand till 1000 nt.
                new_start = int(row[10]) - number_add_length  # New coordinates for Start
                new_end = int(row[11]) + number_add_length  # New coordinates for End

                seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + row[1] + " -range " + str(new_start) + "-" + str(new_end) + " -strand plus -outfmt %s", shell=True, universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
                seq = seq.strip()  # Eliminar EoL caracteres

                new_row = [row[0], row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                chrX_1000nt.append(new_row)

            elif "minus" in row[14]:
                seq_length = int(row[10]) - int(row[11])  # Al reves al ser minus. Podria hacaerlo sino en absoluto, pero bueno, he elegido esta opcion

                number_add_length = int((1000 - seq_length)/2)
                new_start = int(row[10]) + number_add_length  # Al reves al ser minus
                new_end = int(row[11]) - number_add_length  # Al reves al ser minus

                seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + row[1] + " -range " + str(new_end) + "-" + str(new_start) + " -strand minus -outfmt %s", shell=True, universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
                seq = seq.strip()  # Eliminar EoL caracteres

                new_row = [row[0], row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                chrX_1000nt.append(new_row)

    writing_path_input = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_1000nt.csv"
    csv_creator(writing_path_input, chrX_1000nt)

    return (writing_path_input)  # ##Importante para saber el directorio de este archivo creado

# def specific_sequence_1000nt(path_input, chromosome_ID, main_folder_path)

    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    # Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa