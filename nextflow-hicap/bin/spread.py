# spread.py
import pandas as pd
import sys

def spread_data(input_file, output_file):
    df = pd.read_csv(input_file)
    df.columns = df.columns.str.strip()  # Strip leading and trailing spaces from column names
    df_wide = df.pivot(index='SampleFile', columns='Serotype', values='Detection_Status')
    df_wide.reset_index(inplace=True)
    df_wide.to_csv(output_file, index=False)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    spread_data(input_file, output_file)