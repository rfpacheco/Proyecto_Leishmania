import os
import pandas as pd
import subprocess
import csv
import shutil
from pathlib import Path

from modules.aesthetics import boxymcboxface  # Some aesthetics function
from modules.identifiers import genome_specific_chromosome_main
from modules.filters import global_filters_main
from modules.files_manager import folder_creator

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_dic(path_input, path_output):
    """
    Creation af a BLAST database of our whole genome. It uses the BLAST :sup:`R` command line, see BLAST
    `Command Line Application User Manual`_ for more information.


    The generation of the properly database will be placed in the directory where ``path_input`` is.
    It is recommended to use a dedicated folder to this FASTA file so the database is written next to it.

    :param path_input: path to a FASTA file.
    :type path_input: string

    :param path_output: path to the output folder where the BLAST database will be created.
    :type path_output: string

    :return: a BLAST database.
    :rtype: Muitiples files (**.nhr**, **.nin**, **.nog**, **.nsd**, **.nsi** and **.nsq** extensions)
    """

    # Remember is "path.input.dic_path" for "argparse".
    try:
        # boxymcboxface("BLASTn Database creator started")
        cmd = f"makeblastdb -in {path_input} -dbtype nucl -out {path_output}"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)        # print("\nBlast Dictionary created in", path_input)
    except Exception:
        print("\nError: Blast Dictionary couldn't be created")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_blaster(query_path, dict_path, perc_identity):
    """
    This module calls for `blastn` in the BLAST :sup:`R` command line.
    See `Command Line Application User Manual`_ for more information.

    :param query_path: Path to our FASTA query used in BLASTn against our database.
    :type query_path: string

    :param dict_path: Path to the FASTA file ``query_path`` will be launched to. In the same directory should be the BLAST data base created with :func:`~blastn_dic`.
    :type dict_path: string

    :param outfile_path: Path where the results will be saved. Name the whole file, not only the path.
    :type outfile_path: string

    :param perc_identity: Percent of sequence identity which we want to make the BLASTn. **Important**.
    :type perc_identity: string

    :return: A CSV file with the results of the BLASTn.
    :rtype: CSV file


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

    cmd = "blastn -word_size 11 -query " \
        + query_path + " -db " \
        + dict_path \
        + " -perc_identity " + str(perc_identity) \
        + " -outfmt '10 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand sseq'"
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)  # Important the E value
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand", "sseq"]
    return data
    

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def repetitive_blaster(data_input, genome_fasta, folder_path, numbering, maximun_runs):
    """
    This function will iterate till a number of ``maximun_runs`` defined.

    :param genome_fasta: Path to our whole genome sequence in .fasta. It should already be a BLASTn dictionary.
    :type genome_fasta: string

    :param path_input: pandas DataFrame containing data from the first BLASTn. It's the output from :func:`~blastn_blaster`.
    :type path_input: pandas Data Frame

    :param folder_path: Path to a folder where the results will be placed. Subfolder will be created with the cromosome names.
    :type folder_path: string

    :param numbering: Indicates the number showed for the first program run. If it's 0, then, the first run will be name 0.
    :type numbering: integer

    :param maximun_runs: Indicates the last iterative execution of function.
    :type maximun_runs: integer
    """

    # Call the aesthetics function RUN identifier.
    boxymcboxface("RUN " + str(numbering))

    # -----------------------------------------------------------------------------
    # Getting all chromosome IDs from the data_input like "LinJ.01", "LinJ.02", etc. in order.
    chromosome_IDs = sorted(list(data_input.loc[:,"sseqid"].unique()))
    # SHOULD A MAKE A PANDAS GROUPY OBJECT?

    # -----------------------------------------------------------------------------
    for chromosome_ID in chromosome_IDs:
        # if chromosome_ID != "LinJ.01":  # In case we'll need to delete searches of a special chromosome
        genome_specific_chromosome_main(path_input,
                                        chromosome_ID,
                                        folder_path,
                                        genome_fasta)

    # -----------------------------------------------------------------------------
    # Y cuando termine creando el archivo MIXER, lo que hago es purificarlo completamente
    global_filters_main_output = folder_path + "MIXER.csv"  # This one's got the call to "blastn_blaster"
    global_filters_main(global_filters_main_output,
                        global_filters_main_output,
                        genome_fasta)

    # pdb.set_trace()
    folder_output = folder_path + "RUNS"
    folder_creator(folder_output)

    RUN_SAVER_Output = folder_output + "/run_" + str(numbering) + ".csv"
    shutil.copyfile(global_filters_main_output, RUN_SAVER_Output)

    backup_counter = 0
    if backup_counter == 0:
        first_backup = Path(path_input)
        first_backup2 = str(first_backup.parents[0]) + "/" + first_backup.with_suffix("").name + "_backup.csv"
        shutil.copyfile(path_input, first_backup2)  # Here we save the first file just in case
        backup_counter += 1

    shutil.copyfile(global_filters_main_output, path_input)  # This way I restart the first "path_input" with the new "mixed input".
    os.remove(global_filters_main_output)  # Removes "mixer"

    # -----------------------------------------------------------------------------
    if numbering == maximun_runs:
        print("\n\n\nEND of PROGRAM")
    if numbering < maximun_runs:
        numbering += 1
        repetitive_blaster(path_input, genome_fasta, folder_path, numbering, maximun_runs)