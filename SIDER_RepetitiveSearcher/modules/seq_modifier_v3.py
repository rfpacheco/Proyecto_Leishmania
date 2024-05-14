import csv
import subprocess
import re

from modules.files_manager import csv_creator

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def specific_sequence_1000nt(data_input, chromosome_ID, main_folder_path, genome_fasta):
    """
    This function will expand the selected sequences to 1000 nt and will create a CSV file with it. For this, it will use ``blastdbcmd`` from BLAST `Command Line Application User Manual`_, which will get the sequences in the **whole fasta genome** file so every sequence returned will be "plus".


    It uses the function :func:`~modules.file_manager.csv_creator`

    :param path_input: Path to the main CSV file where data will be filtered. In this case is the output of expanded sequences to 1000nt from :func:`~modules.identifiers.specific_sequence_extractor`.
    :type path_input: string

    :param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*, which were present more then once in the ``path_input``. See :func:`~modules.blaster.repetitive_blaster` for more information.
    :type chromosome_ID: member of a python list

    :param main_folder_path: Path where the results will be placed. It will create a subfolder with the chromosome_ID name + "_1000nt.csv".
    :type main_folder_path: string

    :return: a CSV file with the sufix "_1000nt.csv"
    :rtype: CSV file
    """
    # -----------------------------------------------------------------------------
    for index, element in data_input.iterrows():
        if "plus" in element["strand"]:
            lower_coor = element["start"].astype(int)  # We get the start of the sequence
            upper_coor = element["send"].astype(int)  # We get the end of the sequence
        else:  # If it's the "-" strand
            lower_coor = element["send"].astype(int)
            upper_coor = element["start"].astype(int)
        
        subject_length = upper_coor - lower_coor + 1
        if subject_length < 1000:  # If the sequence is less than 1000nt, we'll expand it.
            leftover_length = 1000 - subject_length  # We get the difference between 1000 and the length of the sequence
            leftover_length_halved = leftover_length / 2  # We divide it by 2 to get the half of the difference
            lower_coor = lower_coor - leftover_length_halved  # We subtract the half of the difference to the start
            upper_coor = upper_coor + leftover_length_halved  # We add the half of the difference to the end

            # Check if the coordinates are not out of the genome
            if lower_coor < 0:
                lower_coor = 1
                upper_coor += leftover_length_halved  # Since `new_start` is 1, we add the half of the difference to the end
            
            # We get the sequence from the whole genome
            cmd = f"blastdbcmd -db {genome_fasta} -entry {element['seqid']} -range {lower_coor}-{upper_coor} -strand {element["strand"]} -outfmt %s"
            seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

            # Now with the data let's modify the data frame
            data_input.loc[index, "pident"] = ""
            data_input.loc[index, "length"] = ""
            data_input.loc[index, "qsstart"] = ""
            data_input.loc[index, "qend"] = ""

            if "plus" in element["strand"]:
                data_input.loc[index, "sstart"] = lower_coor
                data_input.loc[index, "send"] = upper_coor
            else:
                data_input.loc[index, "sstart"] = upper_coor
                data_input.loc[index, "send"] = lower_coor

            data_input.loc[index, "evalue"] = ""
            data_input.loc[index, "bitscore"] = ""
            data_input.loc[index, "qlen"] = ""
            data_input.loc[index, "slen"] = len(seq)
            data_input.loc[index, "sseq"] = seq
        else:  # If the sequence is already 1000nt we just pass
            pass
        
    return data_input
            
 

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def specific_sequence_corrected(path_input, nucleotides1000_directory, main_folder_path, chromosome_ID, genome_fasta):
    """
    The main use of this function is to get the real "coordinates" of the sequence.


    1. First we read the 1000nt CSV file. And we extract from row[0] (i.e., Query row) every *different* query (not repeated.
    2. Then we open the 1000nt_BLASTER file from the 1000nt against each other. To understand it:

       - In the 1000nt file there are "x" number of rows. Each row will be given an index like "Seq_z_LinJ.01_plus", were "z" is the index.
       - So the first row (row[0]) in the 1000nt file will be given "Seq_1.." as a name. And every row with the "Seq_1.." subject in the 1000nt_BLASTER file will the result from the first row in the 1000nt file (row[0]).
       - There will be a maximum of "x" "Seq_x.." in the 1000nt_BLASTER file.

    3. In the 1000nt_BLASTER we search every result from a specific "Seq_z..". And we get the maximum alignment and it's coordinates.
    4. We'll get the "z" number and search it's corresponding sequence in the 1000nt file. And with that and the coordinates, we use ``blastcmd`` to get the correct coordinates from the 1000nt sequence. We need to be careful if the sequence is "plus" or "minus".

    :param path_input: Path to the CSV file we'll use to filter data. It's the result from a BLASTn made between the expanded 1000nt sequences.
    :type path_input: string

    :param nucleotides1000_directory: *return* result from the function :func:`~specific_sequence_1000nt`
    :type nucleotides1000_directory: string

    :param main_folder_path: Path where we'll save the CSV data file to.
    :type main_folder_path: string

    :param chromosome_ID: Identification of the cromosome, e.g., "LinJ.07"
    :type chromosome_ID: string

    :return: CSV File with the corrected coordinates of our sequences.
    :rtype: CSV file
    """

    # First from the BLASTn (one againts each other), we get the IDs of the sequences.
    names = []
    with open(path_input, "r") as main_file:
        reader = csv.reader(main_file, delimiter=",")
        for row in reader:
            if row[0] not in names:  # We select the chromosome ID form row[0], so we can have a no-repeated list of those.
                names.append(row[0])  # And example would be "Seq_2_LinJ.01_plus"

    # -----------------------------------------------------------------------------
    # Now, we'll get from the BLASTn (one againts each other), the best alignment.
    # pdb.set_trace()
    chr_x_corrected = []
    for query in names:  # For each chromosome ID row[0] in "names"
        start = []
        end = []
        diference_end_minu_start = 0  # Check it out later. It's to compare it with "difference"

        with open(path_input, "r") as main_file:  # Reads the "_1000nt_Blaster.csv"
            reader = csv.reader(main_file, delimiter=",")
            for row in reader:
                if query in row[1]:  # This is row[1], e.g., "Seq_2_LinJ.01_plus"
                    if row[0] != query:  # This is to remove the sequence that overlaps with itself. So if "eq_2_LinJ.01_plus" overlaps with "eq_2_LinJ.01_plus", we don't analyze it further.
                        # We don't need to difference between "+" and "-" strands, because after BLASTn all the results behave like "+".
                        difference = int(row[11]) - int(row[10])  # Due to how the code is made. Now row[11] will always be > row[10]
                        # Here we'll iterate until we get the larger "difference", i.e., larger "alignment length".
                        if difference > diference_end_minu_start:  # If it's greater than 0
                            diference_end_minu_start = difference  # We save that difference, i.e., alignment length
                            start = []  # Reset "start"
                            end = []  # Reset "end"
                            start.append(int(row[10]))  # Save "start of alignment in subject"
                            end.append(int(row[11]))  # Save "end of alignment in subject"

        # -----------------------------------------------------------------------------
        # Remember we don't need to difference between "+" and "-".
        # In this part, for a specific "query" we'll have the bigger "alignment length" wit its "start" and "end".
        # pdb.set_trace()
        # With the variable "difference" made, this can be removed
        if len(start) > 0 and len(end) > 0:  # I mean this
            min_start = min(start)  # And this
            max_end = max(end)  # And this last one.

            # pdb.set_trace()
            correct_seq = query  # Changed the name to understand it better for the next part.
            number_for_location = re.search("\d+", correct_seq).group()  # Using "regex" to extract the numbers
            number_for_location = int(number_for_location) - 1  # Because Python starts at 0

             #  Now we filter "correct_seq" to obtain a number to filter a CSV list of 4 x 1000 without doing BLAST
            rows_by_number = []  # I need this part to know if it's the correct row while doing comparisons. This way I can compare it with "rows_by_number[0]" or "[4]" or "[3]" without going in order. We do this in the CSV 4 x 1000nt before the blaster to themselves.
            with open(nucleotides1000_directory, "r") as main_file:
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    rows_by_number.append(row)  # Here we get all the rows from the CSV

            # pdb.set_trace()
            with open(nucleotides1000_directory, "r") as main_file:
                reader = csv.reader(main_file, delimiter=",")
                for row in reader:
                    if row == rows_by_number[number_for_location]:  # Asi me aseguro que estoy en la adecuada. Quizas es rizar el rizo pero no se me ocurre en el momento un paso mejor --> ESTO ESTA MAL
                        if "plus" in row[14]:

                            # Since in that file, they are extended to 1000nt, doing 1000 - max end, gives me how much do I have to rest for the sequence end.
                            x = 1000 - max_end
                            new_start = int(row[10]) + min_start - 1
                            new_end = int(row[11]) - x

                            seq = subprocess.check_output("blastdbcmd -db " + genome_fasta + " -entry "
                                                          + row[1] + " -range " + str(new_start) + "-" + str(new_end)
                                                          + " -strand plus -outfmt %s",
                                                          shell=True,
                                                          universal_newlines=True)  # Very important subprocess
                            seq = seq.strip()  # Remove EoL characteres

                            new_row = [query, row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                            chr_x_corrected.append(new_row)

                        elif "minus" in row[14]:  # Remember by how are the coordinates in the "Minus" strand, which are in the position 3'---> 5'
                            x = 1000 - max_end
                            new_start = int(row[10]) - min_start + 1
                            new_end = int(row[11]) + x
                            if new_end <= 0: new_end = 1

                            seq = subprocess.check_output("blastdbcmd -db " + genome_fasta + " -entry "
                                                          + row[1] + " -range " + str(new_end) + "-" + str(new_start)
                                                          + " -strand minus -outfmt %s",
                                                          shell=True,
                                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
                            seq = seq.strip()  # Eliminar EoL caracteres

                            new_row = [query, row[1], "", str(len(seq)), row[4], row[5], "", "", "", "", str(new_start), str(new_end), "", "", row[14], seq]

                            chr_x_corrected.append(new_row)

        # pdb.set_trace()
        if len(start) == 0 and len(end) == 0:  # For the cases where the homology is with itself, we discard it. I may be better to change the code in the future to insert these sequences.
            print("\nALERT: individual " + query + " has no homology with no other seq, so it will not be added to the corrected seqs")

    writing_path_input = main_folder_path + chromosome_ID + "/" + chromosome_ID + "_Corrected.csv"
    csv_creator(writing_path_input, chr_x_corrected)

    return (writing_path_input)