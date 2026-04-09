import sys
import csv

def parse_samples(csv_file):
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        found_header = False
        for row in reader:
            if not row or len(row) < 2:
                continue
            if row[0] == 'File_1':
                found_header = True
                continue
            if not found_header:
                continue
            
            # If we reach here, we are in the data rows
            r1 = row[0]
            r2 = row[1]
            # Generate a sample name by removing _R1_001.fastq.gz from the first file
            sample_name = r1.replace('_R1_001.fastq.gz', '')
            print(f"{sample_name} {r1} {r2}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python parse_samples.py samples.csv")
    else:
        parse_samples(sys.argv[1])
