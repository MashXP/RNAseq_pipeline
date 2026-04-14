import sys

# Usage: python gtf2bed.py input.gtf > output.bed
# A simple script to convert GTF to BED12.

def parse_gtf(gtf_file):
    transcripts = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue
            
            chrom = parts[0]
            start = int(parts[3]) - 1 # GTF is 1-based, BED is 0-based
            end = int(parts[4])
            strand = parts[6]
            
            # Extract transcript_id
            attr_parts = parts[8].split(';')
            transcript_id = None
            for attr in attr_parts:
                if 'transcript_id' in attr:
                    transcript_id = attr.split('"')[1]
                    break
            
            if transcript_id:
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = []
                transcripts[transcript_id].append((chrom, start, end, strand))

    # Write BED12
    for tid, exons in transcripts.items():
        exons.sort() # Sort by start position
        chrom = exons[0][0]
        strand = exons[0][3]
        tx_start = exons[0][1]
        tx_end = exons[-1][2]
        
        block_sizes = []
        block_starts = []
        
        for exon in exons:
            block_sizes.append(str(exon[2] - exon[1]))
            block_starts.append(str(exon[1] - tx_start))
            
        print(f"{chrom}\t{tx_start}\t{tx_end}\t{tid}\t0\t{strand}\t{tx_start}\t{tx_end}\t0\t{len(exons)}\t{','.join(block_sizes)},\t{','.join(block_starts)},")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python gtf2bed.py input.gtf")
    else:
        parse_gtf(sys.argv[1])
