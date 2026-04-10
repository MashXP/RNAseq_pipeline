import pandas as pd
import sys
import os

def parse_featurecounts(input_file, output_file):
    """
    Parses a featureCounts output file (biotype counts) and formats it for MultiQC.
    """
    # Read the featureCounts output, skipping the header comments
    try:
        df = pd.read_csv(input_file, sep="\t", comment="#")
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return
    
    # The first column is the Geneid (or Biotype in our case)
    # The columns from index 6 onwards are the sample counts
    
    if 'Geneid' not in df.columns:
        print(f"Error: 'Geneid' column not found in {input_file}")
        return

    # Set Geneid (Biotype) as index
    df = df.set_index('Geneid')
    
    # Select only the count columns (usually start from column 6, after Length)
    count_cols = df.columns[5:] 
    
    # Create a new DataFrame for MultiQC
    mqc_df = df[count_cols].copy()
    
    # Clean up sample names in columns
    mqc_df.columns = [os.path.basename(c).replace("_Aligned.sortedByCoord.out.bam", "") for c in mqc_df.columns]
    
    # Transpose so samples are rows and biotypes are columns (MultiQC format)
    mqc_df = mqc_df.transpose()
    
    # Save to file with MultiQC custom content header
    with open(output_file, 'w') as f:
        f.write("# id: 'biotype_counts'\n")
        f.write("# section_name: 'Biotype Counts'\n")
        f.write("# description: 'Reads mapped to different gene biotypes (e.g. protein_coding, rRNA).'\n")
        f.write("# plot_type: 'bargraph'\n")
        f.write("# pconfig:\n")
        f.write("#    id: 'biotype_counts_plot'\n")
        f.write("#    title: 'Biotype Distribution'\n")
        f.write("#    ylab: 'Read Counts'\n")
        f.write("#    cpswitch_counts_label: 'Number of Reads'\n")
        f.write(mqc_df.to_csv(sep="\t", index=True, index_label="Sample"))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python biotype_to_multiqc.py input_counts.txt output_mqc.txt")
        sys.exit(1)
        
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    parse_featurecounts(input_path, output_path)
    print(f"Generated MultiQC custom content: {output_path}")
