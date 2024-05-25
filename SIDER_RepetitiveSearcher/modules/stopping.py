import pandas as pd

# Compare the coordinates
def coincidence_counter (df1, df2):
    coincidence = 0
    for _, row in df2.iterrows():
        if row["sstart"] in list(df1["sstart"].values) and row["send"] in list(df1["send"].values):
            coincidence += 1
    return coincidence

def stopping_main(data_df1, data_df2):
    data_df1_plus = data_df1[data_df1["sstrand"] == "plus"].copy()
    data_df1_minus = data_df1[data_df1["sstrand"] == "minus"].copy()
    data_df2_plus = data_df2[data_df2["sstrand"] == "plus"].copy()
    data_df2_minus = data_df2[data_df2["sstrand"] == "minus"].copy()

    coincidence_plus = coincidence_counter(data_df1_plus, data_df2_plus)
    coincidence_minus = coincidence_counter(data_df1_minus, data_df2_minus)

    total_coincidence = coincidence_plus + coincidence_minus
    print(f"\t\t\t- Coincidence with last corrected sequences: {total_coincidence}/{data_df2.shape[0]}")

    if total_coincidence == data_df2.shape[0]:
        print(f"\t\t\t\t TRUE")
        return True
    else:  # If the the coincidence is not the 100%
        print(f"\t\t\t\t FALSE")
        return False