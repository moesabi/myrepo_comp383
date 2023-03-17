# myrepo_comp383 
# HCMV Transcriptome Analysis Pipeline
This is a pipeline for analyzing the transcriptome of Human Cytomegalovirus (HCMV) using RNA sequencing data.

# Step 1: Download the Data
The first step is to download the RNA sequencing data for four samples from the Sequence Read Archive (SRA) using wget.
```python
import os

os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045')
```

# Step 2: Preprocess the Data

The next step is to preprocess the data using fastq-dump to convert the SRA files to FASTQ format.

```python
os.system('fastq-dump -I --split-files SRR5660030')
os.system('fastq-dump -I --split-files SRR5660033')
os.system('fastq-dump -I --split-files SRR5660044')
os.system('fastq-dump -I --split-files SRR5660045')
```

# Step 3: Build the Reference Genome
We will use the HCMV genome as the reference genome for our analysis. We will download the genome from NCBI and use bowtie2-build to build the index.
```python
import os
import Bio
from Bio import Entrez

Entrez.email = "mohammed16alsawi@gmail.com"
handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype='fasta')
fasta = handle.read()
handle.close()
with open('HCMV.fasta', 'w') as f:
    f.write(fasta)

os.system('bowtie2-build HCMV.fasta HCMV')

```

# Step 4: Map the Reads to the Reference Genome
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
            flag_val

```
# Step 5: Assembling the Longest Contig and Running BLAST+ Search



 ```python
 
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Read the SPAdes assembly output and find the longest contig
assembly_file = "contigs.fasta"  # Replace with your SPAdes output file
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
Perform BLAST search on longest contig from SPAdes assembly

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
