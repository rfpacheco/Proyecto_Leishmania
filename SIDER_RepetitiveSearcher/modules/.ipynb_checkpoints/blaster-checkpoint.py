import os


# 1) Creacion del diccionario Blast

# Creation of the Blast Dic of our whole genome
def BlastN_Dic (path_input):
    """
    Creation af a BLAST database of our whole genome. It uses the BLAST\ :sup:`R` \ command line, see BLAST
    `Command Line Application User Manual for more information`_.


    The generation of the properly database will be placed in the **.fasta directory**.
    It is recommended to use a dedicated folder to this genome.

    :param path_input: path to where our file genome sequence .fasta is
    :type path_input: string


    .. _Command Line Application User Manual for more information: https://www.ncbi.nlm.nih.gov/books/NBK279690
    """

    try:
        os.system("makeblastdb -in " + path_input + " -dbtype nucl -parse_seqids")
        print("\nBlast Dictionary created in", path_input)
    except:
        print("\nError: Blast Dictionary couldn't be created")

# BlastN_Dic(Path_Input)

    #Arg 0: STRING. Directorio del archivo FASTA del cual queremos realizar el diccionario BLAST. Los archivos generados se colocaran en ese misma directorio, por ello se recomienda que este dentro de una carpeta unicamente dedicada a estos archivos.