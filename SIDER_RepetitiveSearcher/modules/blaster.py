import os
import pdb


def blastn_dic(path_input):
    """
    Creation af a BLAST database of our whole genome. It uses the BLAST\ :sup:`R` \command line, see BLAST
    `Command Line Application User Manual for more information`_.


    The generation of the properly database will be placed in the **.fasta directory**.
    It is recommended to use a dedicated folder to this genome.

    :param path_input: path to where our file genome sequence .fasta is
    :type path_input: string


    .. _Command Line Application User Manual for more information: https://www.ncbi.nlm.nih.gov/books/NBK279690
    """
    pdb.set_trace()
    print(path_input)
    directory = os.getcwd()
    print(directory)

    os.system("makeblastdb -in " + path_input.dic_path + " -dbtype nucl -parse_seqids")
    # Remember ".dic_path" is the argument in my argparse_main

    try:
        print(path_input)
        os.system("makeblastdb -in " + path_input.dic_path + " -dbtype nucl -parse_seqids")
        print("\nBlast Dictionary created in", path_input)
    except:
        print("\nError: Blast Dictionary couldn't be created")
