import pandas as pd
import os

from modules.bedops import bedops_coincidence
from modules.files_manager import df_columns_restore, columns_to_numeric

def compare_main(df_1, df_2, folder_path, genome_fasta):
    # -----------------------------------------------------------------------------
    # Select data
    df_1_plus = df_1[df_1["sstrand"] == "plus"].copy()
    df_1_minus = df_1[df_1["sstrand"] == "minus"].copy()
    df_2_plus = df_2[df_2["sstrand"] == "plus"].copy()
    df_2_minus = df_2[df_2["sstrand"] == "minus"].copy()
    # -----------------------------------------------------------------------------
    # Prepare paths
    plus_path = os.path.join(folder_path, "bedops_coincidence_plus")
    minus_path = os.path.join(folder_path, "bedops_coincidence_minus")
    # -----------------------------------------------------------------------------
    # Call BEDOPS on plus
    coincidence_plus, new_data_plus = bedops_coincidence(df_1_plus, df_2_plus, plus_path, "plus", genome_fasta)
    if not new_data_plus.empty:  # If the data frame is not empty
        new_data_plus = df_columns_restore(new_data_plus, df_1)  # Restore the columns
        new_data_plus = columns_to_numeric(new_data_plus)  # Convert columns to numeric
    if not coincidence_minus.empty:  # If the data frame is not empty
        coincidence_plus = df_columns_restore(coincidence_plus, df_1)  # Restore the columns
        coincidence_plus = columns_to_numeric(coincidence_plus)  # Convert columns to numeric
    # -----------------------------------------------------------------------------
    # Call BEDOPS on minus. Special case, because BEDOPS reads the coordinates like the "+" strand.
    # First modify the coordinates.
    df_1_minus[["sstart", "send"]] = df_1_minus[["send", "sstart"]]
    df_2_minus[["sstart", "send"]] = df_2_minus[["send", "sstart"]]

    # And now call BEDOPS on minus
    coincidence_minus, new_data_minus = bedops_coincidence(df_1_minus, df_2_minus, minus_path, "minus", genome_fasta)
    # Restore the coordinates
    if not new_data_minus.empty:  # If the data frame is not empty
        new_data_minus[["sstart", "send"]] = new_data_minus[["send", "sstart"]]  # restore "new_data_minus" coordinates
        new_data_minus = df_columns_restore(new_data_minus, df_1)  # Restore the columns
        new_data_minus = columns_to_numeric(new_data_minus)
    if not coincidence_minus.empty:
        coincidence_minus[["sstart", "send"]] = coincidence_minus[["send", "sstart"]]
        coincidence_minus = df_columns_restore(coincidence_minus, df_1)
    # -----------------------------------------------------------------------------
    # Concatenate both Data Frames
    if not new_data_plus.empty and not new_data_minus.empty:
        new_data = pd.concat([new_data_plus, new_data_minus], ignore_index=True)
    else:  # If both Data Frames are empty
        new_data = pd.DataFrame()

    if not coincidence_plus.empty and not coincidence_minus.empty:
        coincidence_data = pd.concat([coincidence_plus, coincidence_minus], ignore_index=True)  # joins both Data Frames
    else:  # If both Data Frames are empty
        coincidence_data = None # return None
    # -----------------------------------------------------------------------------
    return coincidence_data, new_data

    
    


