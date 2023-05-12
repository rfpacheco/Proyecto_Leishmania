import csv
import subprocess

# from modules.filters import chromosome_filter  # Don't call --> ciruclar import
from modules.files_manager import csv_creator
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# 1)En este caso hacemos lo mismo que antes con los solapamientos pero a escala de todo el genoma. El principal cambio importante y MUY importante es el de la edicion del array chromosome_rows, el cual va hasta la primera funcion de "genome_solap_location_filter"


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.1)Aqui metemos todos los START y END para tenerlos localizados.
def genome_solap_location_filter(chromosome_rows):  # Todo STRING menos chromosome_rows
    """
    """
    pos_plus_start = []  # start position for "plus" strand
    pos_plus_end = []
    pos_minus_start = []
    pos_minus_end = []

    for row in chromosome_rows:
        if "plus" in row[14]:  # row[14] is "minus" or "plus". Concretamente esta parte
            pos_plus_start.append(int(row[10]))
            pos_plus_end.append(int(row[11]))
        else:
            pos_minus_start.append(int(row[10]))
            pos_minus_end.append(int(row[11]))

    return (pos_plus_start, pos_plus_end, pos_minus_start, pos_minus_end)

# genome_solap_location_filter(chromosome_rows)
    # Arg 0: Array con los datos a filtrar


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.2)Aqui Metemos los START y END dentro de grupos por proximidad de coordenadas. Usamos la funcion anterior, a la que le impone el chromosome_rows
def genome_solap_location_grouping(chromosome_rows, DNA_sense, max_diff):  # Todo STRING menos Number, max_diff y chromosome_rows
    """
    """
    # Number es: 0 para Plus:Start, 1 para Plus:End, 2 para Minus:Start, 3 para Minus:End
    position_list = genome_solap_location_filter(chromosome_rows)
    if DNA_sense == "plus":
        position_list_main = [position_list[0], position_list[1]]
    elif DNA_sense == "minus":
        position_list_main = [position_list[2], position_list[3]]

    matrix1 = []
    matrix2 = []

    for main_list in position_list_main:
        if main_list == position_list_main[0]:
            matrix = matrix1
        elif main_list == position_list_main[1]:
            matrix = matrix2
        # De esta forma me aseguro que cada main_list vaya a una matrix diferente y al final las junto

        for position in main_list:
            main_statement = False
            for group in matrix:
                for member in group:
                    if abs(member - position) <= max_diff:
                        group.append(position)  # If it's near, we'll append it to the group
                        main_statement = True
                        break  # Romperia codigo y saldriamos del loop "for member" para coontinuar a "if not Found_Matrix"
            if not main_statement:  # Es lo mismo que si "If Found_Matrix == False". Si group esta vacio entonces vendra aqui directamente a poner el numero en ([]), es decir, en un array 3D
                matrix.append([position])  # If not found, We create a List inside a list [[x]], which I called before as "group".

    matrix_main = [matrix1, matrix2]

    return (matrix_main)

# genome_solap_location_grouping(chromosome_rows, DNA_sense, max_diff)
    # Arg 0: Array a filtrar.
    # Arg 1: STRING. Sentido de la hebra, puede ser "plus" o "minus"
    # Arg 2: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.3)Aqui obtenemos los minimos y maximos dependiendo de la hebra que tengamos. Se usa la funcion anterior y le impone el Genome_Rows
def genome_solap_minmax(chromosome_rows, max_diff):  # Todo STRING menos max_diff
    """
    """
    plus = genome_solap_location_grouping(chromosome_rows, "plus", max_diff)
    minus = genome_solap_location_grouping(chromosome_rows, "minus", max_diff)

    plus_min = [min(x) for x in plus[0]]
    plus_max = [max(x) for x in plus[1]]
    minus_max = [max(x) for x in minus[0]]
    minus_min = [min(x) for x in minus[1]]

    return (plus_min, plus_max, minus_max, minus_min)

# genome_solap_minmax(chromosome_rows, max_diff)
    # Arg 0: Array a filtrar.
    # Arg 1: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.4)Modulo para solapantes, es MUY importante, ya que con el se filtran los solapantes parciales que no han sido filtrados por el modulo de solapantes totales de genome_solap_main. Los analiza por parejas.
def genome_solap_by_pairs(rows_to_filter):  # Su argumento es un ARRAY 3D
    """
    """
    rows_final = []
    for first, second in zip(*[iter(rows_to_filter)] * 2):
        two_sequence_rec = []
        two_sequence_rec.append(first)
        two_sequence_rec.append(second)

        sequence_start = []
        sequence_end = []
        homology1 = []
        e_value1 = []
        bit_score1 = []
        for sequence in two_sequence_rec:

            sequence_start.append(sequence[10])
            sequence_end.append(sequence[11])
            homology1.append(float(sequence[2]))  # Esta y las dos de abajo estan puestas por intentar mantener estos datos, pero en realidad no haria nada de falta
            e_value1.append(float(sequence[12]))
            bit_score1.append(float(sequence[13]))

        homology2 = str(round((homology1[0] + homology1[1]) / 2, 3))
        e_value2 = (e_value1[0] + e_value1[1]) / 2
        e_value2 = str("{:.2e}".format(e_value2))
        bit_score2 = str(round((bit_score1[0] + bit_score1[1]) / 2, 1))

        if "plus" in first[14] and abs(int(first[10]) - int(second[10])) <= 1000:  # Este numero es VITAL
            min_start = min(sequence_start)
            max_end = max(sequence_end)
            seq_length = str(int(max_end) - int(min_start) + 1)

            seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry "
                                          + first[1] + " -range " + min_start + "-" + max_end
                                          + " -strand plus -outfmt %s",
                                          shell=True,
                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
            seq = seq.strip()  # Eliminar EoL caracteres

            new_row = [first[0], first[1], homology2, seq_length, first[4], first[5], "", "", "", "", str(min_start), str(max_end), e_value2, bit_score2, first[14], seq]

            rows_final.append(new_row)

        elif "minus" in first[14] and abs(int(first[10]) - int(second[10])) <= 1000:
            max_start = max(sequence_start)
            min_end = min(sequence_end)
            seq_length = str(int(max_start) - int(min_end) + 1)

            seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry "
                                          + first[1] + " -range " + min_end + "-" + max_start
                                          + " -strand minus -outfmt %s",
                                          shell=True,
                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
            seq = seq.strip()  # Eliminar EoL caracteres

            new_row = [first[0], first[1], homology2, seq_length, first[4], first[5], "", "", "", "", str(max_start), str(min_end), e_value2, bit_score2, first[14], seq]

            rows_final.append(new_row)

    return (rows_final)

# genome_solap_by_pairs(rows_to_filter)

    # Arg 0: array que se analizara por parejas obtenido de genome_solap_main


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.5)Este es el principal --> utiliza todos los de arriba y le lleva chromosome_rows
def genome_solap_main(genome_fasta, naming_short, path_input, max_diff, writing_path_input):  # Todo STRING menos max_diff
    """

    :param genome_fasta:
    :type genome_fasta:

    :param naming_short:
    :type naming_short:

    :param path_input:
    :type path_input:

    :param max_diff:
    :type max_diff:

    :param writing_path_input:
    :type writing_path_input:
    """
    genome_solap_main_matrix = []

    chromosome_number = chromosome_filter(genome_fasta, naming_short)

    for chromosome in chromosome_number:

        solap_main_matrix = []

        chromosome_rows = []  # ESENCIAL PARA LA PRIMERA DEFINICION A LA QUE SE LLAMA. Tiene que ir antes de minmax
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:  # chromosome FILTER
                    chromosome_rows.append(row)

        minmax = genome_solap_minmax(chromosome_rows, max_diff)  # Aqui metemos los minimos y maximos anteriores en n array 3D. [0] --> Plus_min | [1] --> Plus_max | [2] --> Minus_max | [3] --> Minus_min
        plus_start_matrix = []
        plus_end_matrix = []
        minus_start_matrix = []
        minus_end_matrix = []

        # En esta parte se busca las secuencias que tengan tanto el minimo como el maximo, de tal forma, se obtiene el mayor alineamiento.
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:  # chromosome FILTER
                    if "plus" in row[14]:
                        if int(row[10]) in minmax[0] and int(row[11]) in minmax[1]:
                            solap_main_matrix.append(row)
                            plus_start_matrix.append(int(row[10]))
                            plus_end_matrix.append(int(row[11]))
                    elif "minus" in row[14]:
                        if int(row[10]) in minmax[2] and int(row[11]) in minmax[3]:
                            solap_main_matrix.append(row)
                            minus_start_matrix.append(int(row[10]))
                            minus_end_matrix.append(int(row[11]))

        # Ahora, por si hay solapantes teniendo uno el minimo y otra secuencia el maximo, necesitaremos las dos. Por eso, revisamos sobre las coordenadas de lo anterior (para no repetir) y, buscamos esos solapantes.
        solap_segments = []  # Aqui meteria todos los segmentos solapados pequeñitos
        solap_segments_plus_start = []
        solap_segments_plus_end = []
        solap_segments_minus_start = []
        solap_segments_minus_end = []
        with open(path_input, "r") as main_file:
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if chromosome in row[1]:
                    if "plus" in row[14]:
                        if int(row[10]) not in plus_start_matrix and int(row[10]) in minmax[0]:
                            if int(row[10]) not in solap_segments_plus_start:
                                solap_segments.append(row)
                                solap_segments_plus_start.append(int(row[10]))

                        if int(row[11]) not in plus_end_matrix and int(row[11]) in minmax[1]:
                            if int(row[11]) not in solap_segments_plus_end:
                                solap_segments.append(row)
                                solap_segments_plus_end.append(int(row[11]))

                    elif "minus" in row[14]:
                        if int(row[10]) not in minus_start_matrix and int(row[10]) in minmax[2]:
                            if int(row[10]) not in solap_segments_minus_start:
                                solap_segments.append(row)
                                solap_segments_minus_start.append(int(row[10]))

                        if int(row[11]) not in minus_end_matrix and int(row[11]) in minmax[3]:
                            if int(row[11]) not in solap_segments_minus_end:
                                solap_segments.append(row)
                                solap_segments_minus_end.append(int(row[11]))

        solap_by_pairs_definitive = genome_solap_by_pairs(solap_segments)  # Esto es CLAVE
        solap_main_matrix += solap_by_pairs_definitive

        genome_solap_main_matrix += solap_main_matrix

    csv_creator(writing_path_input, genome_solap_main_matrix)

# genome_solap_main(genome_fasta, naming_short, path_input, max_diff, writing_path_input)

    # Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    # Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos.
    # Arg 3: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    # Arg 4: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Recordar la extension .csv
