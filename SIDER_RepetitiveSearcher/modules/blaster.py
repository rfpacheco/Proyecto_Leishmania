import pdb # In case of debbuging
import os
import csv
import shutil

from modules.aesthetics import boxymcboxface  # Some aesthetics function

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

    try:
        boxymcboxface("BLASTn Database creator started")
        os.system("makeblastdb -in " + path_input.dic_path + " -dbtype nucl -parse_seqids")
        print("\nBlast Dictionary created in", path_input)
    except Exception:
        print("\nError: Blast Dictionary couldn't be created")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_blaster(query_path, dict_path, outfile_path, perc_identity):
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
        boxymcboxface("BLASTn searcher initiated")
        os.system("blastn -word_size 28 -query "
                  + query_path + " -db "
                  + dict_path + " -out "
                  + outfile_path + " -perc_identity "
                  + perc_identity + " -outfmt '10 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sstrand sseq'")
        print("\nBlaster succesful", outfile_path, "created.")
    except Exception:
        print("\nError: Blaster couldn't be loaded, somthing happened")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def repetitive_blaster(genome_fasta, path_input, folder_path, naming_short, max_diff, numbering, maximun_runs):

    boxymcboxface("RUN ", + str(numbering))

    # -----------------------------------------------------------------------------
    # Searching all "Linj.XX" in the .fasta file
    chr_IDs = []
    with open(genome_fasta, "r") as main_file:
        for line in main_file:
            if ">" in line:
                chr_IDs.append(line[1:8])  # We this we get "LinjJ.XX", being X a number.

    # -----------------------------------------------------------------------------
    # Searching all "Linj.XX" in the .csv file
    csv_IDs_all = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            csv_IDs_all.append(row[1])  # This way we get all "LinJ.XX" in the second column of the .csv. There will be a lot of repeated ones, but with "chr_IDs" we'll know which ones.

    # -----------------------------------------------------------------------------
    # With "chr_in_objetive" I'll get all chr with more than 1 repeated IDs. i.e., at least 2 with the same name.
    chr_in_objetive = []
    for chromosome in chr_IDs:
        counter = 0
        for objetive in csv_IDs_all:
            if chromosome in objetive:
                counter += 1
        if counter > 1:  # Here we'll add all chrs with more than 1 representative, i.e., repeated
            chr_in_objetive.append(chromosome)

    # -----------------------------------------------------------------------------
    for chromosome_ID in chr_in_objetive:
        # if chromosome_ID != "LinJ.01":  # In case we'll need to delete searches of a special chromosome
        Genome_Specific_chromosome_Main(path_input,
                                        chromosome_ID,
                                        folder_path,
                                        genome_fasta,
                                        naming_short,
                                        max_diff)

    # -----------------------------------------------------------------------------
    # Y cuando termine creando el archivo MIXER, lo que hago es purificarlo completamente
    Global_Filters_Main_Output = folder_path + "/MIXER.csv"
    Global_Filters_Main(Global_Filters_Main_Output,
                        Global_Filters_Main_Output,
                        genome_fasta,
                        naming_short,
                        max_diff)

    RUN_SAVER_Output = folder_path + "/RUNS/run_" + str(numbering) + ".csv"
    shutil.copyfile(Global_Filters_Main_Output, RUN_SAVER_Output)
    shutil.copyfile(Global_Filters_Main_Output, path_input)  # ## Asi reseteo el Path input con el nuevo documento para lanzarlo todo de nuevo. El path input antiguo ya no existe porque ha sido sobre escrito (aunque se ha guardado en RUNS)
    os.remove(Global_Filters_Main_Output)  # ##Eliminamos el Mixer, para que luego se cree de nuevo

    # -----------------------------------------------------------------------------
    if numbering == maximun_runs:
        print("\n\n\nEND of PROGRAM")
    if numbering < maximun_runs:
        numbering += 1
        repetitive_blaster(genome_fasta, path_input, folder_path, naming_short, max_diff, numbering, maximun_runs)