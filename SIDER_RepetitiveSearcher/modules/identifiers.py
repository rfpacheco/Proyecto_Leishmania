import os
import csv

from modules.files_manager import folder_creator, csv_creator


def specific_sequence_extractor(path_input, chromosome_ID, main_folder_path):
    """
    It reads
    """
    chr_x_seqs = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if chromosome_ID in row[1]:
                chr_x_seqs.append(row)

    folder_path = main_folder_path + "/" + chromosome_ID
    folder_creator(folder_path)

    writing_path_input = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + ".csv"

    csv_creator(writing_path_input, chr_x_seqs)

    return (folder_path, writing_path_input)  # Es importante porque asi nos devuelven los nuevos directorios. No se las puedo añadir a variables globales, porque lamentablemente Python no funciona asi

# specific_sequence_extractor(path_input, chromosome_ID, main_folder_path)
    # Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    # Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    # Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa. Lo utilizara para generar una subcarpeta con el nombre del cromosoma


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def genome_specific_chromosome_main(path_input, chromosome_ID, main_folder_path, genome_fasta, naming_short, max_diff):

    New_Directories = specific_sequence_extractor(path_input, chromosome_ID, main_folder_path)
    folder_path = New_Directories[0]  # El directorio de la carpeta del cromosoma
    Last_Output = New_Directories[1]  # El directorio del comando anterior. Estos dos pasos los hago por facilitar la lectura del codigo

    Nucleotides1000_Directory = Specific_Sequence_1000nt(Last_Output, chromosome_ID, main_folder_path)

    Fasta_Creator_Output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_1000nt.fasta"
    Fasta_Creator(Nucleotides1000_Directory, Fasta_Creator_Output)

    BlastN_Dic(Fasta_Creator_Output)

    Blaster_Output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_1000nt_Blaster.csv"
    BLASTN_Blaster(Fasta_Creator_Output,
                   Fasta_Creator_Output,
                   Blaster_Output,
                   "85")

    Filter_by_Column(Blaster_Output,
                     "length",
                     100,
                     Blaster_Output)

    Corrected_Sequences = Specific_Sequence_Corrected(Blaster_Output, Nucleotides1000_Directory, main_folder_path, chromosome_ID)

    Subfamilies_File_Path_Writing = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_Subfamily.csv"
    Subfamily_Sorter(Blaster_Output, Corrected_Sequences, Subfamilies_File_Path_Writing)

    Second_Fasta_Creator_Output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_Corrected.fasta"
    Fasta_Creator(Corrected_Sequences, Second_Fasta_Creator_Output)

    Second_Blaster_Output = main_folder_path + "/" + chromosome_ID + "/" + chromosome_ID + "_BLAST_MAIN.csv"
    BLASTN_Blaster(Second_Fasta_Creator_Output,
                   genome_fasta,
                   Second_Blaster_Output,
                   "60")

    Global_Filters_Main(Second_Blaster_Output,
                        Second_Blaster_Output,
                        genome_fasta,
                        naming_short,
                        max_diff)


    CSV_Mixer_Output = main_folder_path + "/" + "MIXER.csv"
    if os.path.isfile(CSV_Mixer_Output) is False:  # #Cuando no existe, se crea
        CSV_Mixer(path_input, Second_Blaster_Output, CSV_Mixer_Output)  # Para mezclar
    else:  # Si existe ya el archivo porque ha sido creado, se cambia el path_input por CSV_Mixer_Output
        CSV_Mixer(CSV_Mixer_Output, Second_Blaster_Output, CSV_Mixer_Output)

#genome_specific_chromosome_main(path_input, chromosome_ID, main_folder_path, genome_fasta, naming_short, max_diff)

    #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa
    #Arg 4: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    #Arg 5: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 6: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE