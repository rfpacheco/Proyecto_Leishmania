from Bio import SeqIO


# This code let's me choose which sequences I want in a fasta output (NEED TO CHANGE FOR MY USE)
# This one chooses from "../ref/L_infantum_ALL_36Chr.fasta" the chromosome 2 and 3 and outputs it like "../Prueba1/Results/Test1"
# from Bio import SeqIO
with open("../Prueba1/Results/Test1", "w") as out_file:
    for count, fasta in enumerate(SeqIO.parse(open("../ref/L_infantum_ALL_36Chr.fasta"), "fasta"), 1):
        name, sequence = fasta.id, str(fasta.seq)
        if count in (2, 3):
            SeqIO.write(fasta, out_file, "fasta")