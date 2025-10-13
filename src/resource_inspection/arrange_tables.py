import pandas as pd


def safe_to_timedelta(time_str):
    if '-' in time_str:
        days_part, time_part = time_str.split('-')
        return pd.to_timedelta(f"{days_part} days {time_part}")
    parts = time_str.split(':')
    if len(parts) == 2:
        return pd.to_timedelta(f"0:{time_str}")
    return pd.to_timedelta(time_str)


def process_format_1(file_path: str, data_type: str, algorithm: str) -> pd.DataFrame:
    """
    Processes TSV files according to specified rules, generalized for multiple algorithms.
    """
    df = pd.read_csv(file_path, sep='\t')

    # Process IND rows
    ind_df = df[df['job name'] == f"{algorithm}_IND"].copy()
    ind_df['tissue'] = ind_df['job id'].apply(
        lambda x: 'brain' if any(x.endswith(f"_{i}.batch") for i in range(1, 6)) else 'kidney'
    )

    # Process CONCAT and CONCAT_QUANT rows
    concat_df = df[df['job name'].isin([f"{algorithm}_CONCAT", f"{algorithm}_CONCAT_QUANT"])].copy()

    def assign_concat_tissue(row, concat_df):
        if len(concat_df[concat_df['job name'] == row['job name']]) == 2:
            return 'brain' if row['job id'].endswith("_1.batch") else 'kidney'
        elif len(concat_df[concat_df['job name'] == row['job name']]) > 2: # for FLAIR where J&C parallelizes by chromosome
            concat_tissues = sorted(set(concat_df['job id'].str.extract(r'(\d+)_')[0].astype(int)))
            min_tissue = min(concat_tissues)
            return 'brain' if row['job id'].startswith(f"{min_tissue}_") else 'kidney'
        else:
            raise ValueError("Not enough CONCAT jobs to determine tissues.")

    concat_df['tissue'] = concat_df.apply(assign_concat_tissue, axis=1, concat_df=concat_df)

    # Process TAMA_CONDITION rows
    tama_condition_df = df[df['job name'] == "TAMA_CONDITION"].copy()
    tama_condition_df['tissue'] = tama_condition_df['job id'].apply(
        lambda x: 'brain' if x.endswith("_1.batch") else 'kidney'
    )

    # Convert elapsed and totalcpu to timedelta safely
    for df_subset in [ind_df, concat_df, tama_condition_df]:
        df_subset['elapsed'] = df_subset['elapsed'].apply(safe_to_timedelta)
        df_subset['totalcpu'] = df_subset['totalcpu'].apply(safe_to_timedelta)

    # Aggregate IND and TAMA_CONDITION together
    combined_individual_df = pd.concat([ind_df, tama_condition_df])
    aggregated_ind_df = combined_individual_df.groupby('tissue').agg({
        'elapsed': 'sum',
        'totalcpu': 'sum',
        'reqcpus': 'first',
        'maxrss': 'max'
    }).reset_index()
    aggregated_ind_df['strategy'] = "Call&Join"
    aggregated_ind_df['entry'] = "summary"

    # Aggregate CONCAT and CONCAT_QUANT together
    aggregated_concat_df = concat_df.groupby('tissue').agg({
        'elapsed': 'sum',
        'totalcpu': 'sum',
        'reqcpus': 'first',
        'maxrss': 'max'
    }).reset_index()
    aggregated_concat_df['strategy'] = "Join&Call"
    aggregated_concat_df['entry'] = "summary"

    # Add common columns to aggregated dataframes
    for aggregated_df in [aggregated_ind_df, aggregated_concat_df]:
        aggregated_df['data_type'] = data_type
        aggregated_df['algorithm'] = algorithm

    # Set individual entry types
    ind_df['entry'] = "individual"
    ind_df['strategy'] = "Call&Join"
    ind_df['data_type'] = data_type
    ind_df['algorithm'] = algorithm

    tama_condition_df['entry'] = "TAMA"
    tama_condition_df['strategy'] = "Call&Join"
    tama_condition_df['data_type'] = data_type
    tama_condition_df['algorithm'] = algorithm

    concat_df['entry'] = concat_df['job name'].apply(
        lambda x: "quant" if x == f"{algorithm}_CONCAT_QUANT" else "individual"
    )
    concat_df['strategy'] = "Join&Call"
    concat_df['data_type'] = data_type
    concat_df['algorithm'] = algorithm

    # Select relevant columns
    columns_to_keep = ['data_type', 'tissue', 'algorithm', 'strategy', 'entry', 'elapsed', 'reqcpus', 'maxrss', 'totalcpu']
    final_df = pd.concat([
        ind_df[columns_to_keep],
        tama_condition_df[columns_to_keep],
        aggregated_ind_df[columns_to_keep],
        concat_df[columns_to_keep],
        aggregated_concat_df[columns_to_keep]
    ], ignore_index=True)

    return final_df


def process_format_isoseq(file_path: str, data_type: str, algorithm: str) -> pd.DataFrame:
    """
    Processes TSV files according to specified rules for isoseq format.
    """
    df = pd.read_csv(file_path, sep='\t')

    ind_pattern = f"^{algorithm}.*_IND$"
    ind_df = df[df['job name'].str.contains(ind_pattern, regex=True)].copy()
    # Process IND rows
    ind_df = df[df['job name'] == f"{algorithm}_IND"].copy()
    ind_df['tissue'] = ind_df['job id'].apply(
        lambda x: 'brain' if any(x.endswith(f"_{i}.batch") for i in range(1, 6)) else 'kidney'
    )

    # Process CONCAT rows
    concat_pattern = f"^{algorithm}.*_CONCAT$"
    concat_df = df[df['job name'].str.contains(concat_pattern, regex=True)].copy()
    concat_df['tissue'] = concat_df['job id'].apply(
        lambda x: 'brain' if x.endswith("_1.batch") else 'kidney'
    )

    # Process TAMA_CONDITION rows
    tama_condition_df = df[df['job name'] == "TAMA_CONDITION"].copy()
    tama_condition_df['tissue'] = tama_condition_df['job id'].apply(
        lambda x: 'brain' if x.endswith("_1.batch") else 'kidney'
    )

    # Convert elapsed and totalcpu to timedelta safely
    for df_subset in [ind_df, concat_df, tama_condition_df]:
        df_subset['elapsed'] = df_subset['elapsed'].apply(safe_to_timedelta)
        df_subset['totalcpu'] = df_subset['totalcpu'].apply(safe_to_timedelta)

    # Aggregate IND and TAMA_CONDITION together
    combined_individual_df = pd.concat([ind_df, tama_condition_df])
    aggregated_ind_df = combined_individual_df.groupby('tissue').agg({
        'elapsed': 'sum',
        'totalcpu': 'sum',
        'reqcpus': 'first',
        'maxrss': 'max'
    }).reset_index()
    aggregated_ind_df['strategy'] = "Call&Join"
    aggregated_ind_df['entry'] = "summary"

    # Aggregate CONCAT and CONCAT_QUANT together
    aggregated_concat_df = concat_df.groupby('tissue').agg({
        'elapsed': 'sum',
        'totalcpu': 'sum',
        'reqcpus': 'first',
        'maxrss': 'max'
    }).reset_index()
    aggregated_concat_df['strategy'] = "Join&Call"
    aggregated_concat_df['entry'] = "summary"

    # Add common columns to aggregated dataframes
    for aggregated_df in [aggregated_ind_df, aggregated_concat_df]:
        aggregated_df['data_type'] = data_type
        aggregated_df['algorithm'] = algorithm

    # Set individual entry types
    ind_df['entry'] = "individual"
    ind_df['strategy'] = "Call&Join"
    ind_df['data_type'] = data_type
    ind_df['algorithm'] = algorithm

    tama_condition_df['entry'] = "TAMA"
    tama_condition_df['strategy'] = "Call&Join"
    tama_condition_df['data_type'] = data_type
    tama_condition_df['algorithm'] = algorithm

    concat_df['entry'] = concat_df['job name'].apply(
        lambda x: "sqanti" if "SQANTI" in x else ("filter" if "filter" in x else "individual")
    )
    concat_df['strategy'] = "Join&Call"
    concat_df['data_type'] = data_type
    concat_df['algorithm'] = algorithm

    # Select relevant columns
    columns_to_keep = ['data_type', 'tissue', 'algorithm', 'strategy', 'entry', 'elapsed', 'reqcpus', 'maxrss', 'totalcpu']
    final_df = pd.concat([
        ind_df[columns_to_keep],
        tama_condition_df[columns_to_keep],
        aggregated_ind_df[columns_to_keep],
        concat_df[columns_to_keep],
        aggregated_concat_df[columns_to_keep]
    ], ignore_index=True)

    return final_df



def assemble_dataframe(job_details: dict[str, dict[str, str]]) -> pd.DataFrame:
    """
    Processes algorithm entries and returns a combined DataFrame.
    """
    df_list = []

    for data_type, algorithms in job_details.items():
        for algorithm, file_path in algorithms.items():
            print(f"Processing {data_type}: {algorithm} ...")
            if algorithm == "ISOSEQ":
                df = process_format_isoseq(file_path, data_type, algorithm)
            else:
                df = process_format_1(file_path, data_type, algorithm)
            df_list.append(df)

    combined_df = pd.concat(df_list, ignore_index=True)
    return combined_df


if __name__ == "__main__":
    # Example usage:
    job_details: dict[str, dict[str, str]] = {
        "ONT": {
            "ISOQUANT": "/home/fabianje/repos/documenting_NIH/fabian/reports/ont/isoquant/run2_report/job_details.tsv",
            "FLAIR": "/home/fabianje/repos/documenting_NIH/fabian/reports/ont/flair_ar_sr/run2_report/job_details.tsv",
            "bambu": "/home/fabianje/repos/documenting_NIH/fabian/reports/ont/bambu/run2_report/job_details.tsv",
            "TALON": "/home/fabianje/repos/documenting_NIH/fabian/reports/ont/talon/run1_report/job_details.tsv"
        },
        "IsoSeq": {
            "ISOQUANT": "/home/fabianje/repos/documenting_NIH/fabian/reports/isoseq/isoquant/run3_report/job_details.tsv",
            "FLAIR": "/home/fabianje/repos/documenting_NIH/fabian/reports/isoseq/flair_ar_sr/run7_report/job_details.tsv",
            "bambu": "/home/fabianje/repos/documenting_NIH/fabian/reports/isoseq/bambu/run3_report/job_details.tsv",
            "TALON": "/home/fabianje/repos/documenting_NIH/fabian/reports/isoseq/talon/run3_report/job_details.tsv",
            "MANDALORION": "/home/fabianje/repos/documenting_NIH/fabian/src/resource_inspection/alejandro_jobs/mandalorion/job_details.tsv",
            "ISOSEQ": "/home/fabianje/repos/documenting_NIH/fabian/src/resource_inspection/alejandro_jobs/isoseq/job_details.tsv"
        }
    }

    df = assemble_dataframe(job_details)
    # Print all rows and columns
    df.to_csv("test.csv", index=False)
    df[df["entry"] == "summary"].to_csv("summary.csv", index=False)
