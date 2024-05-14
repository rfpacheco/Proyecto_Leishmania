import argparse
import os
import shutil
import time  # to measure the time of the program
import subprocess  # to call the command line

from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.aesthetics import boxymcboxface
from modules.files_manager import fasta_creator

# Initiate parser
parser = argparse.ArgumentParser(
    prog="SIDER_RepetitiveSearcher",
    description="This is a program to search for repetitive sequences in SIDERs elements in Leishmania spp.",
)

# Let's get the user input data
parser.add_argument("-d", "--data", type=str, required=True, help="Path to the input data file")
parser.add_argument("-g", "--genome", type=str, required=True, help="Path to the genome file")

# Parsing the arguments
args = parser.parse_args()

# =============================================================================
# Creating working folder place 
# =============================================================================
# Ask the user for the folder name and data location
folder_name = input("Enter folder name: ")
data_location = input("Enter path where you want to place all data: ")

data_location = os.path.normpath(data_location)  # Normalize the path to avoid problems with the OS

folder_location = os.path.join(data_location, folder_name)  # Create the folder location

# Create the folder with the given name
os.makedirs(folder_location, exist_ok=True)
print(f"Folder {folder_name} created in {data_location}")

# =============================================================================
# Take original data used and copy it inside the folder
# =============================================================================
# Create subdirectory for original data
original_data_folder = os.path.join(folder_location, "original_data")
os.makedirs(original_data_folder, exist_ok=True)

# Copy the data and genome files to the original_data folder
shutil.copy(args.data, original_data_folder)
shutil.copy(args.genome, original_data_folder)
# print(f"Files copied to {original_data_folder}")
args_data_path = os.path.join(original_data_folder, os.path.basename(args.data))  # save the path so we can use this one instead of the original one
args_genome_path = os.path.join(original_data_folder, os.path.basename(args.genome))  # save the path so we can use this one instead of the original one

# =============================================================================
# First blaster automatization
# =============================================================================
# Create folder for main BLASTN dictionary
blastn_dict_path = os.path.join(folder_location, "dict_data")
os.makedirs(blastn_dict_path, exist_ok=True)
blastn_dict_path_out = os.path.join(blastn_dict_path, os.path.basename(args.genome))

# Create the BLASTn dictionary
blastn_dic(args.genome, blastn_dict_path_out)

# Call the first BLASTn
identity_1 = input("Enter the identity for the first BLASTn step: ")  # user input for the identity

boxymcboxface("First BLASTn step initiated")

tic = time.perf_counter()  # Start the timer
first_blaster = blastn_blaster(args_data_path, blastn_dict_path_out, identity_1)  # It has the data frame for the first blaster
toc = time.perf_counter()  # Stop the timer

print(f"==>First BLASTn step took {toc - tic:0.2f} seconds")
print(f"==>First BLASTn row length: {first_blaster.shape[0]}")


# =============================================================================
# Call the second and last BLASTn
# =============================================================================
# 1) Create fasta file for the second BLASTn from first_blaster data frame
fasta_file_path = os.path.join(folder_location, "First_Blaster.fasta")  # Path to the fasta file to create

# Now let's create the fasta file
tic = time.perf_counter()  # Start the timer
fasta_creator(data_input = first_blaster,
              fasta_output_path = fasta_file_path)
toc = time.perf_counter()  # Stop the timer
print(f"==>Fasta file creation took {toc - tic:0.2f} seconds")

repetitive_blaster(data_input = first_blaster,
                   genome_fasta = blastn_dict_path_out,  # path to the genome dict
                   folder_path = folder_location,
                   numbering = 1,
                   maximun_runs = 2)
