import argparse
import os
import shutil
import time  # to measure the time of the program
from datetime import datetime
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
folder_name = input("\n\nEnter folder name: ")
data_location = input("Enter path where you want to place all data: ")

data_location = os.path.normpath(data_location)  # Normalize the path to avoid problems with the OS
folder_location = os.path.join(data_location, folder_name)  # Create the folder location

# Create the folder with the given name
os.makedirs(folder_location, exist_ok=True)
print(f"{"."*20} Folder {folder_name} created in {data_location}")

identity_1 = input("Enter the identity for the first BLASTn step: ")  # user input for the identity

# =============================================================================
# Start time
# =============================================================================
start_time = datetime.now()
tic_main = time.perf_counter()  # Start the timer
formatted_start_time = start_time.strftime("%Y %B %d at %H:%M")
print(f"{"."*20} Program started: {formatted_start_time}")
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
blastn_dic(path_input=args.genome, 
           path_output=blastn_dict_path_out)

# Call the first BLASTn
boxymcboxface(message="First BLASTn step initiated")

tic = time.perf_counter()  # Start the timer
first_blaster = blastn_blaster(query_path=args_data_path, 
                               dict_path=blastn_dict_path_out, 
                               perc_identity=identity_1)  # It has the data frame for the first blaster
first_blaster.to_csv(os.path.join(folder_location, "First_Blaster.csv"), index=False, header=0, sep=",")  # Save the data frame to a CSV file
toc = time.perf_counter()  # Stop the timer
print(f"1. Initial data:\n",
      f"\t- Data row length: {first_blaster.shape[0]}\n",
      f"\t- Execution time: {toc - tic:0.2f} seconds")




# =============================================================================
# Call the second and last BLASTn
# =============================================================================
# 1) Create fasta file for the second BLASTn from first_blaster data frame
fasta_file_path = os.path.join(folder_location, "First_Blaster.fasta")  # Path to the fasta file to create

# Now let's create the fasta file
tic = time.perf_counter()  # Start the timer
fasta_creator(data_input=first_blaster,
              fasta_output_path=fasta_file_path)
toc = time.perf_counter()  # Stop the timer
print("")
print(f"2. Fasta file creation:\n",
      f"\t- Execution time: {toc - tic:0.2f} seconds")

tic = time.perf_counter()  # Start the timer
repetitive_blaster(data_input=first_blaster,
                   genome_fasta=blastn_dict_path_out,  # path to the genome dict
                   folder_path=folder_location,
                   numbering=1,
                   maximun_runs=2,
                   start_time=formatted_start_time)
toc = time.perf_counter()  # Stop the timer

# =============================================================================
# End time
# =============================================================================
toc_main = time.perf_counter()  # Stop the timer
end_time = datetime.now()
formatted_end_time = end_time.strftime("%Y %B %d at %H:%M")
boxymcboxface(message="END OF THE PROGRAM")
print(f"\t- Execution time: {toc_main - tic_main:0.2f} seconds\n",
      f"\t- Program started: {formatted_start_time}\n",
      f"\t- Program ended: {formatted_end_time}")
