# myrepo_comp383 
# HCMV Transcriptome Analysis Pipeline
This is a pipeline for analyzing the transcriptome of Human Cytomegalovirus (HCMV) using RNA sequencing data.

# Download the Data
The first step is to download the RNA sequencing data for four samples from the Sequence Read Archive (SRA) using wget. We are going to use python scripts
```python
import os
# Download SRR5660030..033..44..45 file from the SRA database
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030') #convertsfrompythontobeUNIX
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045')
```

# Preprocess the Data

The next step is to preprocess the data using fastq-dump to convert the SRA files to FASTQ format.

```python

import os 

#split the separate forward and reverse reads and convert the specified SRA files into FASTQ format using -I and --split


os.system('fastq-dump -I --split-files SRR5660030')
os.system('fastq-dump -I --split-files SRR5660033')
os.system('fastq-dump -I --split-files SRR5660044')
os.system('fastq-dump -I --split-files SRR5660045')
```

# Build the Reference Genome
We will use the HCMV genome as the reference genome for our analysis. We will download the genome from NCBI and use bowtie2-build to build the index.
```python
import os
import Bio
from Bio import Entrez
# Fetch the HCMV genome sequence from the NCBI nucleotide database
Entrez.email = "mohammed16alsawi@gmail.com"
handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype='fasta')
fasta = handle.read()
handle.close()
with open('HCMV.fasta', 'w') as f:
    f.write(fasta)
# Build the Bowtie2 index using the HCMV genome sequence
os.system('bowtie2-build HCMV.fasta HCMV')
# Run Bowtie2 alignment for each of the four samples
print('starting bowtie2 .sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam')

```

# Map the Reads to the Reference Genome
We will use bowtie2 to map the reads to the reference genome and filter out the unmapped reads.
```python
import os

# Getting initial read counts
initial_counts = {}
samples_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
for s in samples_list:
    with open(f'{s}_1.fastq') as file1:
        initial_counts[s] = sum([1 for line in file1]) / 4
    print(f'{s} has {initial_counts[s]:,} read pairs before filtering')

# Running Bowtie2 on samples and keeping only HCMV index mapped reads
for s in samples_list:
    output_sam = f'{s}_mapped.sam'
    os.system(f'bowtie2 --quiet -x HCMV -1 {s}_1.fastq -2 {s}_2.fastq -S {output_sam}')

    # Filtering unmapped reads by reading the SAM file and writing mapped reads to a FASTQ file
    post_filter_counts = 0
    with open(output_sam) as samfile, open(f'{s}_HCMV.fastq', 'w') as fastqfile:
        for line in samfile:
            if line.startswith('@'):
                continue

            parts = line.split('\t')
            flag_val = int(parts[1])

            if flag_val & 4 == 0:
                post_filter_counts += 1
                qname_val = parts[0]
                seq_val = parts[9]
                qual_val = parts[10]
                fastqfile.write(f"@{qname_val}\n{seq_val}\n+\n{qual_val}\n")

    post_filter_counts /= 4
    print(f'{s} has {post_filter_counts:,} read pairs after filtering')

    # Writing read counts to a file
    with open('log.txt', 'a') as logfile:
        logfile.write(f'{s} has {initial_counts[s]:,} read pairs before filtering and {post_filter_counts:,} read pairs after filtering.\n')






```
# SPAdes assembly 
We  use SPAdes to perform a genome assembly of four different transcriptomes (SRR5660030, SRR5660033, SRR5660044, SRR5660045) and also append the SPAdes command to a log file
```python
import os

samples = {
    'SRR5660030': ('SRR5660030_HCMV.fastq'),
    'SRR5660033': ('SRR5660033_HCMV.fastq'),
    'SRR5660044': ('SRR5660044_HCMV.fastq'),
    'SRR5660045': ('SRR5660045_HCMV.fastq'),
}

# Assemble all four transcriptomes together using SPAdes
spades_input = ''
for index, (sample_name, fq) in enumerate(samples.items(), start=1):
    spades_input += f'--s{index} {fq} '

spades_command = f'spades.py -k 77,99,127 -t 4 --only-assembler {spades_input.strip()} -o HCMV_SRR_assembly'
os.system(spades_command)

# Write the SPAdes command to the log file
with open('log.txt', 'a') as log_file:
    log_file.write(f'SPAdes command: {spades_command}\n')
    
    
```



# Assembling the Longest Contig and Running BLAST+ Search



 ```python
 
 from Bio import SeqIO

# Input file
assembly_file = 'contigs.fasta'

# Initialize counters
contigs_gt_1000 = 0
assembly_length = 0

# Iterate through contigs in the assembly file
for record in SeqIO.parse(assembly_file, 'fasta'):
    contig_length = len(record.seq)
    
    # Check if contig length is greater than 1000
    if contig_length > 1000:
        contigs_gt_1000 += 1
        assembly_length += contig_length

# Write results to the log file
with open('log.txt', 'a') as log_file:
    log_file.write(f'There are {contigs_gt_1000} contigs > 1000 bp in the assembly.\n')
    log_file.write(f'There are {assembly_length} bp in the assembly.\n')

```
# Analyzing HCMV Transcriptomes Using SPAdes and BLAST
reads the SPAdes assembly output and identifies the longest contig which is then saved to a file. The longest contig is then used in a BLAST+ NCBI nr  
```python
 
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Read the SPAdes assembly output and find the longest contig
assembly_file = "contigs.fasta"  
contigs = SeqIO.parse(assembly_file, "fasta")
longest_contig = max(contigs, key=lambda x: len(x))

# Save the longest contig to a file
with open("longest_contig.fasta", "w") as output_handle:
    SeqIO.write(longest_contig, output_handle, "fasta")

# Use the longest contig as a BLAST+ input to query the nr nucleotide database
query_sequence = str(longest_contig.seq)
blast_program = "blastn"
blast_database = "nr"

search_parameters = {
    "entrez_query": "txid10357[Organism:exp]",  # Betaherpesvirinae taxonomy ID
    "hitlist_size": 10,  # Number of results to return
}

result_handle = NCBIWWW.qblast(
    blast_program,
    blast_database,
    query_sequence,
    **search_parameters
)

# Save the result to a file
with open("blast_result.xml", "w") as output_file:
    output_file.write(result_handle.read())

result_handle.close()
```

# Perform BLAST search on longest contig from SPAdes assembly
outputs the top 10 outputs with the criteria, "sacc   pident   length   qstart   qend   sstart   send   bitscore   evalue   stitle"

```python

from Bio.Blast import NCBIXML

# Parse the BLAST result XML file
with open("blast_result.xml", "r") as result_handle:
    blast_records = list(NCBIXML.parse(result_handle))

# Write the top 10 hits to a log file
with open("blast_top_hits.log", "w") as log_file:
    # Write the header row
    log_file.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

    # Iterate over the BLAST records and alignments
    for record in blast_records:
        for alignment in record.alignments[:10]:  # Limit to top 10 hits
            hsp = alignment.hsps[0]  # Only keep the best HSP

            # Extract the information and format the output
            output_line = f"{alignment.accession}\t{hsp.identities * 100 / hsp.align_length:.2f}\t{hsp.align_length}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{hsp.bits}\t{hsp.expect}\t{alignment.title}\n"

            # Write the output line to the log file
            log_file.write(output_line)

            ```
