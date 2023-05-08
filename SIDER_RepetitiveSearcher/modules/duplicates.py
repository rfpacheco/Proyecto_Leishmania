import csv

from modules.filters import chromosome_IDs
from modules.files_manager import csv_creator


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1) Este filtro de duplicado ya funcionan a nivel del genoma entero. Se revisan las secuencias de las filas una,y por el cromosoma, sentido de la hebra y agurpaciones, se eliminan las duplicaciones.

# 1.1) El primero es el filtro para Duplicados, funciona igual que el anterior pero es para todo un genoma y necesita tanto el CSV como el fasta.


def genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense, max_diff):  # Todo STRING menos max_diff

    matrix_all_genome = []
    chromosome_number = chromosome_IDs(genome_fasta, naming_short)  # I obtain a Python list, e.g., ["LinJ.01", "LinJ.02", ...]

    for chromosome in chromosome_number:  # Aqui nos metemos dentro del cromosoma a buscar.
        location_start = []  # Almacenaremos los START de ese cromosoma
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")  # Here we read the CSV file "BLAST_MAIN.csv"
            for row in reader:
                if chromosome in row[1]:  # CHROMOSOME FILTER
                    if DNA_sense in row[14]:  # DNA_sense hay que poner PLUS o MINUS
                        location_start.append(int(row[10]))  # We save the "Start of alignment in query" from the CSV file.

            # -----------------------------------------------------------------------------
            matrix_filter = []
            position_global = []  # #!!!!!MUY IMPORTANTE¡¡¡¡Con esto le decimos que no repita localizaciones, es vital
            with open(path_input, "r") as main_file:  # we read the CSV "BLAST_MAIN.csv" again
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    if chromosome in row[1]:  # CHROMOSOME FILTER
                        if DNA_sense in row[14] and int(row[10]) not in position_global:
                            # Arriba le hemos dicho que revise las localizaciones en el "position_global", para que no se repita. Se hace aquí porque "position_rec" va cambiando constantemente.
                            position_rec = []
                            for position in location_start:
                                if abs(int(row[10]) - position) < max_diff:  # In this part we make sure we are NEAR our position in the genome. If the "position" is 35000 and row[10] is 35400 and max_diff = 600, then abs(35400-35000)=400 < 600, so we are near.
                                    # If position is 35000 and row[10] is 10000, then abs(10000-35000)=25000 > 600, so we are FAR AWAY.
                                    if position not in position_rec:  # Esta parte es para que en lo que dura una organizacion de "position_rec", no se repitan ningun valor dentro de ella y por consiguiente, dentro del Global
                                        position_rec.append(position)
                                        position_global.append(position)

                            DNAseq_filter = []
                            with open(path_input, "r") as main_file:  # Tengo que abrirlo de nuevo para empezar en la primera row siempre
                                reader = csv.reader(main_file, delimiter=",")
                                for row in reader:
                                    if chromosome in row[1]:  # CHROMOSOME FILTER
                                        if int(row[10]) in position_rec:  # Comprobamos si la localizacion de START esta dentro de "position_rec" de este momento -recordar que va cambiando-. Si esta dentro de estas localizaciones, entonces nos centramos en mirar duplicaciones solo dentro de estas localizaciones.
                                            if row[15] in DNAseq_filter:
                                                continue  # Si la secuencia esta dentro de nuestra base de datos, el CONTINUE salta el resto del codigo y vuelve al anterior loop FOR, evitando asi añadir duplicados
                                            else:  # Al no estar la secuencia dentro de nuestra base de datos, se le añade a nuestra base de datos Y a la matrix final global.
                                                DNAseq_filter.append(row[15])
                                                matrix_filter.append(row)
        matrix_all_genome += matrix_filter

    return (matrix_all_genome)

# genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense, max_diff)

    # Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el FASTA del genoma entero
    # Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 3: STRING. Sentido de la hebra, puede ser "plus" o "minus"
    # Arg 4: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE


# 1.2) Este es el principal:


def genome_duplicate_filter(genome_fasta, naming_short, path_input, max_diff, writing_path_input):  # Todo STRING menos max_diff

    DNA_sense = ["plus", "minus"]

    matrix_main1 = genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense[0], max_diff)

    matrix_main2 = genome_pre_duplicate_filter(genome_fasta, naming_short, path_input, DNA_sense[1], max_diff)

    matrix_main = matrix_main1 + matrix_main2

    csv_creator(writing_path_input, matrix_main)

# genome_duplicate_filter(genome_fasta, naming_short, path_input, max_diff, writing_path_input)

    # Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    # Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 3: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    # Arg 4: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, Recordar la extension .csv