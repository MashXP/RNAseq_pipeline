import sys
import csv

def parse_samples(csv_file):
    count = 0
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        found_header = False
        for row in reader:
            # Skip pre-header rows
            if not found_header:
                if row and row[0].strip() == 'File_1':
                    found_header = True
                continue

            # Stop at ANY empty or non-FASTQ row after the header
            # (handles the "Genome files:" section and trailing blank lines)
            if not row or not row[0].strip() or not row[0].strip().endswith('.fastq.gz'):
                break

            r1 = row[0].strip()
            r2 = row[1].strip()
            sample_name = r1.replace('_R1_001.fastq.gz', '')
            print(f"{sample_name} {r1} {r2}")
            count += 1

    if count == 0:
        print("ERROR: parse_samples.py found 0 samples. Check CSV format and 'File_1' header.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python parse_samples.py samples.csv")
        sys.exit(1)
    else:
        parse_samples(sys.argv[1])
