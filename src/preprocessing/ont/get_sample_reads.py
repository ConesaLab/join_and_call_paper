import os
import pandas as pd
import glob

def extract_lengths(file_path):
    try:
        df = pd.read_csv(file_path, delimiter='\t') 
        fourth_column = df["length"].tolist()
        return(fourth_column)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

def main(directory_path):
    file_data = {}
    for file_path in glob.glob(directory_path):
        column_data = extract_lengths(file_path)
        if column_data is not None:
            file_data[file_path.split("/")[3]] = column_data


    df = pd.DataFrame({ key:pd.Series(value) for key, value in file_data.items() })
    #df.to_csv("../analysis/read_lengths_df.csv")
    df.to_csv("../analysis/read_lengths_SIRV_df.csv")

if __name__ == '__main__':
    #directory_path = "../analysis/run_SQANTI/*/*_classification.txt"
    directory_path = "../analysis/run_SQANTI/*_SIRV/*_classification.txt"
    
    main(directory_path)
