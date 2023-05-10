# At the moment this module it's due to be modified
import csv

from modules.files_manager import csv_creator

# 4)Este apartado simplemente nos generara un archivo con las subfamilias para ver como se clasifican


def subfamily_sorter(path_input, corrected_elements_path, writing_path_input):
    """
    .. warning::
       Due to be modified
    """

    names = []
    with open(path_input, "r") as main_file:  # It opens the 1000nt BALSTn results after "filter_by_column"
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if row[0] not in names:
                names.append(row[0])

    for_subfamilies = []
    with open(path_input, "r") as main_file:  # ##path_input seria el documento del Blaster
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            new_row = [row[0], row[1], row[2], "", ""]  # importante las posiciones "" vacias, porque ahi añadiremos las longitudes para hacer el calculo del siguiente filtrado por longitud
            for_subfamilies.append(new_row)

    for_subfamilies2 = []
    with open(corrected_elements_path, "r") as main_file:  # ##Aqui tendriamos el documento de los elementos ya corregidos.
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            new_row = [row[0], row[3]]
            for_subfamilies2.append(new_row)

    for row_sub2 in for_subfamilies2:
        for row_sub1 in for_subfamilies:
            if row_sub2[0] in row_sub1[0]:
                row_sub1[3] = row_sub2[1]

            if row_sub2[0] in row_sub1[1]:
                row_sub1[4] = row_sub2[1]

    # Aqui vamos a cambiar el % a valores de 1.
    for_subfamilies3 = []
    for row in for_subfamilies:
        if row[3] != "" and row[4] != "":  # Para eliminar los que no hacen homologia con ninguno y tienen por tanto los valores vacios
            subfamily_filter = ((float(row[2]) * int(row[3])) / int(row[4])) / 100

            if 0.85 <= subfamily_filter <= 1.15:
                new_row = [row[0], row[1]]
                for_subfamilies3.append(new_row)

        if row[3] == "" and row[4] == "":  # #Para el que no hace homologia con ninguno salvo consigo mismo, al menos que en el archivo de subffamilias, aparezca, pero solo.
            subfamily_filter = 1.0
            new_row = [row[0], row[1]]
            for_subfamilies3.append(new_row)

    for_subfamilies4 = []
    rec_for_subfamilies = []
    for sequence in names:
        for pair in for_subfamilies3:
            if sequence in pair[0]:
                if pair[1] not in rec_for_subfamilies:
                    rec_for_subfamilies.append(pair[1])

        for_subfamilies4.append(sorted(rec_for_subfamilies))
        rec_for_subfamilies = []

    # Ahora nos encargaremos de eliminar los duplicados:
    for_subfamilies5 = []
    for groups in for_subfamilies4:
        if groups not in for_subfamilies5:
            for_subfamilies5.append(groups)

    csv_creator(writing_path_input, for_subfamilies5)

# subfamily_sorter(path_input, corrected_elements_path, writing_path_input)
    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: RETURN del resultado de la funcion Specific_sequence_Corrected
    # Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa