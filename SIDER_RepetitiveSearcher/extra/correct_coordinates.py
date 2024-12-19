import argparse  # type: ignore

import pandas as pd  # type: ignore
import subprocess  # type: ignore
import os  # type: ignore
from Bio import SeqIO  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore

from SIDER_RepetitiveSearcher.extra.functions import blastn_dic, blastn_blaster


# noinspection DuplicatedCode
def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter elements based on chromosome count and e-value.")
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('-d', '--dict_path', type=str, required=True, help='Path to the genome fasta file.')
    parser.add_argument('--word_size', type=str, required=True, help='word_size value of BLASTN')
    parser.add_argument('--recaught_file', type=str, required=True, help='Path to the recaught file.')
    return parser.parse_args()

def true_sider_test(data, genome_dict, evalue):
    matches = pd.Series([False] * data.shape[0])
    not_matches = pd.Series([False] * data.shape[0])
    accepted = 0
    rejected = 0
    for index, row in data.iterrows():
        print("="*50)
        print(f"Analyzing row {index + 1} of {data.shape[0]}")

        # Create fasta temp file:
        fasta_id = f"Seq_{index}"
        fasta_seq = row["sseq"]
        fasta_query = f"<(echo -e '>{fasta_id}\n{fasta_seq}')"

        # BLASTn the fasta tmp file against the genome
        blastn_data = blastn_blaster(fasta_query, genome_dict, evalue, args.word_size)

        # Test SIDER
        # noinspection DuplicatedCode
        if blastn_data.count("\n") <= 1:
            not_matches[index] = True
            rejected += 1
            print("\t\tREJECTED")
        else:
            blastn_data = blastn_data.strip().split("\n")
            blast_data_df = pd.DataFrame([x.split(",") for x in blastn_data if x])
            if blast_data_df[1].nunique() >= 5:
                matches[index] = True
                accepted += 1
                print("\t\tACCEPTED")
            else:
                not_matches[index] = True
                rejected += 1
                print("\t\tREJECTED")
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")

    # noinspection DuplicatedCode
    print(f"The total number of matches is: {matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / data.shape[0] * 100, 2)}%")
    print("~"*50)
    print(f"The total number of not matches is: {not_matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of not matches is: {round(not_matches.sum() / data.shape[0] * 100, 2)}%")

    accepted_elem = data[matches]
    rejected_elem = data[~matches]
    return accepted_elem, rejected_elem

if __name__ == "__main__":
    args = parse_arguments()
    file_path = os.path.expanduser(args.file)
    file_parent = os.path.dirname(file_path)
    folder_path = os.path.join(file_parent, 'correct_coordinates')
    os.makedirs(folder_path, exist_ok=True)

    # read the data
    df = pd.read_csv(args.file, sep=",", header=0)

    # Create BLASTn dictionary
    genome_path = os.path.expanduser(args.dict_path)
    genome_folder_path = os.path.join(folder_path, 'blastn_dict')
    os.makedirs(genome_folder_path, exist_ok=True)
    genome_file_path_out = os.path.join(genome_folder_path, os.path.basename(genome_path))
    blastn_dic(genome_path, genome_file_path_out)

    # Apply True SIDER Test
    ## expected value = 1.0E-09
    real_siders, false_siders = true_sider_test(df, genome_file_path_out, 1.0E-09)

