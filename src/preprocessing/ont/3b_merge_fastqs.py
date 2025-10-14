import os

def generate_concatenation_string(base_dir, runs, barcodes):
    concatenation_strings = []

    for barcode in barcodes:
        file_paths = []
        
        for run in runs:
            run_dir = os.path.join(base_dir, f"../data/fastqs/run{run}/fastq_pass/merged_fastqs/")
            
            resqued_input_file = os.path.join(run_dir, f"barcode{barcode:02d}_merged_rescued.fastq")
            full_length_input_file = os.path.join(run_dir, f"barcode{barcode:02d}_merged_full_length.fastq")
            
            if os.path.exists(resqued_input_file):
                file_paths.append(resqued_input_file)
                
            if os.path.exists(full_length_input_file):
                file_paths.append(full_length_input_file)
        
        concatenation_string = " ".join(file_paths)
        concatenation_strings.append(concatenation_string)
        
        print(f"zcat {concatenation_string} > concatenated_barcode{barcode:02d}.fastq")

# Define the base directory
base_directory = '.'  # or the path to the directory containing run2 to run8

# Define the range of runs and barcodes
runs = range(2, 9)  # run2 to run8
barcodes = range(1, 25)  # barcode01 to barcode24

# Generate the concatenation strings
generate_concatenation_string(base_directory, runs, barcodes)