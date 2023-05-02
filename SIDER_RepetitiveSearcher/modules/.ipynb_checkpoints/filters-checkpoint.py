import csv

from modules.file_manager import csv_creator


# 2)Version mejorada para utilizarlo segun la columna que se desee editar, sea por longitud de las secuencias como por porcentaje de homologia

def filter_by_column(path_input, column, size_filter, writing_path_input):  # Todo STRING menos size_filter
    if column == "length":
        column = 3
    elif column == "percent":
        column = 2

    matrix_filter_by_column = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")  # Recordar que antes al poner outfmt 10, ahora estan separados por comas.
        for row in reader:
            if column == 3:
                if 1000 >= int(row[column]) >= size_filter:  # ##1000 por el problema de duplicaciones que me dijo Requena
                    matrix_filter_by_column.append(row)
            elif column == 2:
                if float(row[column]) >= size_filter:  # Necesario para pasar de STRING a FLOAT
                    matrix_filter_by_column.append(row)
    csv_creator(writing_path_input, matrix_filter_by_column)

#filter_by_column(path_input, column, size_filter, writing_path_input)

    #Arg 0: STRING. Directorio del archivo en formato CSV al que queremos filtrar los datos.
    #Arg 1: STRING. Puede ser "length" o "percent", dependiendo de lo que quereamos filtrar.
    #Arg 2: INT. Numero que nos indica el minimo para filtrar dependiendo del Arg 1.
    #Arg 3: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, recordar poner la extension .csv