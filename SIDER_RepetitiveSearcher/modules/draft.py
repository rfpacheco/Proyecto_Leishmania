import csv
import subprocess

from modules.files_manager import csv_creator


def specific_sequence_corrected(path_input, nucleotides1000_directory, main_folder_path, chromosome_ID):
    """
    """

    names = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if row[0] not in names:  # We select the chromosome ID form row[0], so we can have a no-repeated list of those.
                names.append(row[0])

    # -----------------------------------------------------------------------------
    chr_x_corrected = []
    for query in names:  # For each chromosome ID row[0] in "names"
        start = []
        end = []

        diference_end_minu_start = 0  # #Si resulta ser mayor, se le añadira el valor mayor
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if query in row[1]:  # Mejor hacerlo con este para realizar lo de las tablas del CSV que hice
                    if row[0] != query:  # Lo  row[0] != query es para eliminar la secuencia que solapa consigo misma. Por otro lado, no hace falta diferenciar entre plus y minus. Porque hemos lanzado nuestras pequeñas secuencias a un blast propio, y en si se comportan como si todas fueran "plus"
                        difference = int(row[11]) - int(row[10])  # Por como esta el codigo ahora 11 siempre sera mayor que 10
                        if difference > diference_end_minu_start:

                            diference_end_minu_start = difference
                            start = []
                            end = []
                            start.append(int(row[10]))
                            end.append(int(row[11]))

        # De nuevo, no hace falta diferenciar entre plus y minus por lo mismo de antes

        # Esta parte ya no se necesita despues de haber indicado lo de difference, pero en otro momento la quito
        if len(start) > 0 and len(end) > 0:  # Asi creo que evito el hecho de que seq5 no tenga homologia con nada.
            min_start = min(start)
            max_end = max(end)

            correct_seq = query

            mumber_for_location = int(correct_seq[4]) - 1  # Asi cogemos con INT la cuarta posicion de seq_X, siendo X un numero y le restamos 1, porque ya sabemos que en Python todo empieza en 0
            # Ahora filtramos el correct_seq para obtener un numero para poder filtrar la lista de csv de 4 x 1000 sin hacer blaster.

            rows_by_number = []  # Esta parte la necesito para poder saber luego si estoy en la row correcta al hacer comparaciones, ya que puedo compararlo con rows_by_number[0] o [4] o [3], sin ir en orden. Esto se hace en el csv de 4 x 1000nt antes del blaster a ellos mismos
            with open(nucleotides1000_directory, "r") as main_file:
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    rows_by_number.append(row)

            with open(nucleotides1000_directory, "r") as main_file:
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    if row == rows_by_number[mumber_for_location]:  # Asi me aseguro que estoy en la adecuada. Quizas es rizar el rizo pero no se me ocurre en el momento un paso mejor
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