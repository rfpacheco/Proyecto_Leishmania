import pdb
from Bio import SeqIO
# from pathlib import Path


# This thing here to now my "functions" imported or made that doesn't start with "_"
print("\n".join([x for x in dir() if not x.startswith("_")]))


# This code let's me choose which sequences I want in a fasta output (NEED TO CHANGE FOR MY USE)
# This one chooses from "../ref/L_infantum_ALL_36Chr.fasta" the chromosome 2 and 3 and outputs it like "../Prueba1/Results/Test1"
# from Bio import SeqIO
def fasta_extractor(pathfile, outfile, extract_tuple):
    """
    This function let me chose the sequences I want from the fasta by their number.

    Packages needed:
       - ``from Bio import SeqIO``
       - ``from pathlib import Path``

    :param pathfile: Path to the fasta file we want to subset sequences
    :type pathfile: string

    :param outfile: Path to the output file. Include name and ".fasta" suffix.

    :param extract_tuple: typle with our index numbers to extrac the sequence. For example (2, 4, 6) will extract sequences 2, 4 and 6 from ``pathfile``
    :type extract_tuple: tuple

    :return: A fasta file with the sequence we want
    :rtype: fasta file
    """
    with open(outfile, "w") as out_file:
        # Remember "enumerate" starts in "1"
        for count, fasta in enumerate(SeqIO.parse(open(pathfile), "fasta"), 1):  # from Bio import SeqIO
            # name, sequence = fasta.id, str(fasta.seq)
            if count in extract_tuple:
                SeqIO.write(fasta, out_file, "fasta")