import os
import csv
import pdb  # For debugging


def folder_creator(options):
    """
    Creates a Folder in the current path

    :paran options: Name of the folder
    :type options: string

    :returns: Folder in the current directory
    """
    pdb.set_trace()  # Debugging mark
    options = str(options.file_name)
    if not os.path.exists(options):
        os.mkdir(options)
        print("\nDirectory", options, "created at:\n",
              os.path.abspath(options))
    else:
        print("\nDirectory", options, "already exists")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def csv_creator(writing_path_input, writing_input):
    """
    This function will create .csv files

    :param writing_path_input: Path where the .csv file will be created
    :type writing_path_input: string

    :param writing_input: What info we want to add to the .csv file, normally a Matrix/Array 3D
    :type writing_input: TO BE DONE
    """
    with open(writing_path_input, "w") as OutCSV:
        writer = csv.writer(OutCSV)
        writer.writerows(writing_input)
        print("\nCSV:", writing_path_input, "has been created.")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def csv_mixer(path_input1, path_input2, writing_path_input):
    """
    Mixes two .csv files. Then, uses :func:`~csv_creator` to generate a .csv with the mixed files.

    :param path_input1: Path to the first .csv to mix. It will take the first rows in the final mixed file.
    :type path_input1: string

    :param path_input2: Path to the second .csv to mix. It will take rows adther the first .csv.
    :type path_input2: string

    :param writing_path_input: path to both .csv files
    :type writing_path_input: string
    """

    csv_mixer_matrix = []
    with open(path_input1, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            csv_mixer_matrix.append(row)
    with open(path_input2, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            csv_mixer_matrix.append(row)

    csv_creator(writing_path_input, csv_mixer_matrix)
    
