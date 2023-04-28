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

    :param query_path: Path to our .fasta query used in BLASTn against our database.
    :type query_path: string

    :param dict_path: Path to our database we made with :func:`~blastn_dic`.
    :type dict_path: string

    :param outfile_path: path where the results will be saved. **Remember to write it with the .csv ending**.
    :type outfile_path: string

    :param perc_identity: Percent oh identity which we want to make the BLASTn
    :type perc_identity: string
    
    :return: a .csv file with the results of the BLASTN
    :rtype: .csv file


    .. list-table:: Arguments used in ``blastn``
       :header-rows: 1

       * - Argument
         - Type
         - Description

       * - ``-word_size``
         - integer
         - | 28 by default.
           | BLAST search starts with finding a perfect sequence match of length given by this parameter.
           | See `world_size info`_ for more information.

       * - ``-query``
         - string
         - Query file name.

       * - ``-db``
         - string
         - BLAST database name.

       * - ``-out``
         - string
         - Output file name.

       * - ``-perc_identity``
         - integer
         - Percent identity cutoff.

       * - ``-outfmt``
         - string
         - | Alignment view options. Written inside in a nested string.
           | Explained in :ref:`Table <outfmt_Table>`
           | See `Command Line Application User Manual`_ and `Metagenomics BLASTn manual`_ for help.


    ..
        # This is to mark the next table

    .. _outfmt_Table:


    .. list-table:: Arguments used in ``--outfmt`` inside ``blastn``
       :header-rows: 1


       * - Argument
         - Description

       * - ``10``
         - Comma-separated values.

       * - ``qseqid``
         - | Query Seq-id.
           | ID of our query sequence.

       * - ``sseqid``
         - | Subject Seq-id.
           | ID of our subject sequence.

       * - ``pident``
         - Percentage of identical matches.

       * - ``length``
         - Alignment length.

       * - ``qlen``
         - Query sequence length.

       * - ``slen``
         - Subject sequence length.

       * - ``mismatch``
         - Number of mismatches.

       * - ``gapopen``
         - Number of gap openings.

       * - ``qstart``
         - Start of alignment in query.

       * - ``qend``
         - End of alignment in query.

       * - ``sstart``
         - Start of alignment in subject.

       * - ``send``
         - End of alignment in subject.

       * - ``evalue``
         - Expect value.

       * - ``bitscore``
         - Bit score.

       * - ``sstrand``
         - Subject strand.

       * - ``sseq``
         - Aligned part of subject sequence.
    """

    try:
        os.system("blastn -word_size 28 -query "
                  + query_path + " -db "
                  + dict_path + " -out "
                  + outfile_path + " -perc_identity "
                  + perc_identity + " -outfmt '10 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sstrand sseq'")
        print("\nBlaster succesful", outfile_path, "created.")
    except Exception:
        print("\nError: Blaster couldn't be loaded, somthing happened")

