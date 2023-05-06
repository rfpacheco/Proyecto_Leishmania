import csv

from modules.files_manager import csv_creator


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def filter_by_column(path_input, column, size_filter, writing_path_input):
    """
    This function will filter a CSV data depending on ``length`` (if we want to firlter by sequence length) or ``percent`` (if we want to filter by identity percent).

    :param path_input: Path to the .csv file we want to filter data.
    :type path_input: string

    :param column: Can be ``length`` (if we want to firlter by sequence length) or ``percent`` (if we want to filter by identity percent)
    :type column: string

    :param size_filter: Number to filter dependint of the **column** argument.
    :type size_filter: integer

    :param writing_path_input: Path to the CSV file this function will create and save
    :type writing_path_input: string

    :return: A CSV file with the dalta filtered depending on the **column** and **size_filter** argumetns.
    :rtype: CSV file
    """

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

