import os
import pdb

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_dic(path_input):
    """
    Creation af a BLAST database of our whole genome. It uses the BLAST\ :sup:`R` \command line, see BLAST
    `Command Line Application User Manual`_ for more information.


    The generation of the properly database will be placed in the **.fasta directory**.
    It is recommended to use a dedicated folder to this genome.

    :param path_input: path to where our file genome sequence .fasta is
    :type path_input: string


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
    except Exception:
        print("\nError: Blast Dictionary couldn't be created")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blast_blaster(query_path, dict_path, outfile_path, perc_identity):
    """
    This module calls for `blastn` in the BLAST\ :sup:`R` \command line.
    See `Command Line Application User Manual`_ for more information.

    .. list-table:: Arguments used in ``blastn``
       :header-rows: 1

       * - Argument
         - Description
       * - ``word_size``
         - | 11 by default.
           | BLAST search starts with finding a perfect sequence match of length given by this parameter.
           | See `world_size info`_ for more information.
        
        
           
    

    """
    try:
        os.system("blastn -word_size 11 -query "
                  + query_path + " -db "
                  + dict_path + " -out "
                  + outfile_path + " -perc_identity "
                  + perc_identity + " -outfmt '10 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sstrand sseq'")
        print("\nBlaster succesful", outfile_path, "created.")
    except Exception:
        print("\nError: Blaster couldn't be loaded, somthing happened")

