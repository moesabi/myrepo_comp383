# myrepo_comp383 
# HCMV Transcriptome Analysis Pipeline
This is a pipeline for analyzing the transcriptome of Human Cytomegalovirus (HCMV) using specific sequencing data.

# Download the Data
The first step is to download the RNA sequencing data for four samples from the Sequence Read Archive (SRA) using wget. We are going to use python scripts
```python
import os #executes system commands
#download SRR5660030, SRR5660033, SRR5660044, SRR5660045 files from the SRA database
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030') #os. converts from python
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044')
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045')
```

# Preprocess the Data

The next step is to preprocess the data using fastq-dump to convert the SRA files to FASTQ format.

```python

import os 

#split the forward and reverse reads and convert the SRA files into FASTQ format using "-I" and "--split"


os.system('fastq-dump -I --split-files SRR5660030')
os.system('fastq-dump -I --split-files SRR5660033')
os.system('fastq-dump -I --split-files SRR5660044')
os.system('fastq-dump -I --split-files SRR5660045')
```

# Build the Reference Genome
Use the HCMV genome as the reference genome for analysis. Download the genome from NCBI and use bowtie2-build to build the index.
```python
import os  
import Bio
from Bio import Entrez
#fetches HCMV genome sequence from NCBI database
Entrez.email = "mohammed16alsawi@gmail.com" #email for Entrez module
handle = Entrez.efetch(db="nucleotide", id='NC_006273.2', rettype='fasta') #retrieves genome sequence from NCBI 
fasta = handle.read() 
handle.close()
with open('HCMV.fasta', 'w') as f: #save the HCMV genome sequence as a fasta file
    f.write(fasta)
#builds the Bowtie2 index using the HCMV genome sequence
os.system('bowtie2-build HCMV.fasta HCMV')
#run Bowtie2 alignment for all four samples
print('starting bowtie2 .sam') #displays the start of alignment process 'start bowtie2.sam' 
os.system('bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam')

#Output file in GitHub as "log1_reads_transcriptome_b4_Bowtie.txt"

```

# Map the Reads to the Reference Genome
Use bowtie2 to map the reads to the reference genome and filter out the unmapped reads.
```python
import os

# created dictionary to store initial read counts
initial_counts = {}
#creates list of sample names
samples_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

#loop through each sample in the list
for s in samples_list:
    with open(f'{s}_1.fastq') as file1: # opens the file with the sample name and '_1.fastq' ext.
        initial_counts[s] = sum([1 for line in file1]) / 4 #gets the initial read count by counting lines in the file and dividing by 4
    print(f'{s} has {initial_counts[s]:,} read pairs before filtering')

# Running Bowtie2 on samples and keeping only HCMV index mapped reads
for s in samples_list:
    output_sam = f'{s}_mapped.sam'
    os.system(f'bowtie2 --quiet -x HCMV -1 {s}_1.fastq -2 {s}_2.fastq -S {output_sam}')

    # Filtering unmapped reads by reading the SAM file and writing mapped reads to a FASTQ file
    post_filter_counts = 0
    
    #opens the output SAM file and a new FASTQ file for writing
    with open(output_sam) as samfile, open(f'{s}_HCMV.fastq', 'w') as fastqfile: # loops through each line in the SAM file
        for line in samfile:
            if line.startswith('@'): #skips the header lines
                continue

            parts = line.split('\t') #splits lines into columns
            flag_val = int(parts[1]) # gets flag value from the second column

            if flag_val & 4 == 0: # check if the read is mapped ('flag value & 4 is 0')
                post_filter_counts += 1 #runs increments of the count of filtered reads
                qname_val = parts[0] # retrieves the values from the SAM file columns
                seq_val = parts[9]
                qual_val = parts[10]
                fastqfile.write(f"@{qname_val}\n{seq_val}\n+\n{qual_val}\n")  # write the mapped read to the new FASTQ file

    post_filter_counts /= 4 # divides the post-filter count by 4 to get the number of all read pairs
    print(f'{s} has {post_filter_counts:,} read pairs after filtering')

    # appends the read counts to a file
    with open('log.txt', 'a') as logfile: # appends the read counts to a log file
        logfile.write(f'{s} has {initial_counts[s]:,} read pairs before filtering and {post_filter_counts:,} read pairs after filtering.\n')






```
# SPAdes assembly 
Use SPAdes to perform a genome assembly of four different transcriptomes (SRR5660030, SRR5660033, SRR5660044, SRR5660045) and also append the SPAdes command to a log file
```python
import os
# define the dictionary to store sample names and their corresponding FASTQ files
samples = {
    'SRR5660030': ('SRR5660030_HCMV.fastq'),
    'SRR5660033': ('SRR5660033_HCMV.fastq'),
    'SRR5660044': ('SRR5660044_HCMV.fastq'),
    'SRR5660045': ('SRR5660045_HCMV.fastq'),
}

#assemble all four transcriptomes together using SPAdes

spades_input = '' #initializes the empty string to store input arguments for the SPAdes command
for index, (sample_name, fq) in enumerate(samples.items(), start=1): # loop through all of the samples by extracting the sample name and FASTQ file
    # appends the current sample's input flag and FASTQ file to the spades_input string
    spades_input += f'--s{index} {fq} '

spades_command = f'spades.py -k 77,99,127 -t 4 --only-assembler {spades_input.strip()} -o HCMV_SRR_assembly'
os.system(spades_command)

# wite the SPAdes command to the log file
with open('log.txt', 'a') as log_file:
    log_file.write(f'SPAdes command: {spades_command}\n')
    
    
```



# Assembling the Longest Contig and Running BLAST+ Search
Parse the file containing contigs in FASTA format and count the number of contigs longer than 1000 bp and then calculate the total assembly length and write the results to a log file.

 ```python
 
 from Bio import SeqIO

# define the assembly file and makes the input file containing contigs in FASTA format
assembly_file = 'contigs.fasta'

# initialize counters
contigs_gt_1000 = 0 # initialize counter for the number of contigs > 1000 bp
assembly_length = 0 # initialize a variable to store the total assembly length

#Iterate's through contigs in the assembly file
for record in SeqIO.parse(assembly_file, 'fasta'): 
    # calculates the length of the current contig
    contig_length = len(record.seq)
    
    # This checks to see if contig length is greater than 1000 bp
    if contig_length > 1000:
        contigs_gt_1000 += 1  
        assembly_length += contig_length #adds the contig length to the total assembly length

# opens the log file in append mode to write the results
with open('log.txt', 'a') as log_file:
    log_file.write(f'There are {contigs_gt_1000} contigs > 1000 bp in the assembly.\n') #write's the number of contigs longer than 1000 bp
    log_file.write(f'There are {assembly_length} bp in the assembly.\n') # writes the total assembly length


#Output file in GitHub as "log2_nmbr_contigs_nmbr_assmbly.txt"
```
# Analyzing HCMV Transcriptomes Using SPAdes and BLAST

Read the SPAdes assembly output file, find the longest contig, and save to a new file. Then uses the longest contig as a query to perform a BLAST+ search against the nr nucleotide database and saves the results in an XML file.

```python
 
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML #from the BLAST+ ppt. 

#reads the SPAdes assembly output and find the longest contig
assembly_file = "contigs.fasta"  #inputs the assembly file
contigs = SeqIO.parse(assembly_file, "fasta") #parses the input file as FASTA format
longest_contig = max(contigs, key=lambda x: len(x)) #finds the longest contig

#saves the longest contig to a file
with open("longest_contig.fasta", "w") as output_handle: #opens the output file to be appended 
    SeqIO.write(longest_contig, output_handle, "fasta") #writes the longest contig to the output file

#use the longest contig as a BLAST+ input to query the nr nucleotide database
query_sequence = str(longest_contig.seq) # convert the longest contig sequence to a string
blast_program = "blastn" #uses blastn program to use
blast_database = "nr" #searches through nr database

#defines the search parameters
search_parameters = {
    "entrez_query": "txid10357[Organism:exp]",  # this is the betaherpesvirinae taxonomy ID
    "hitlist_size": 10,  # the number of results to return
}

result_handle = NCBIWWW.qblast( #performs BLAST search using the parameters
    blast_program,
    blast_database,
    query_sequence,
    **search_parameters
)

# Save the result to a file
with open("blast_result.xml", "w") as output_file: 
    output_file.write(result_handle.read()) #write the BLAST result to the output file

result_handle.close() #close the result handle 
```

# Perform BLAST search on longest contig from SPAdes assembly

Parse the BLAST result XML file using NCBIXML module and extract the top 10 alignments for each. Write into to a log file, including the accession number, percentage identity, alignment length, query start and end, subject start and end, bitscore, e-value, and alignment title

```python

from Bio.Blast import NCBIXML #biopython library modules

#Opens and parses the BLAST result XML file
with open("blast_result.xml", "r") as result_handle:
    blast_records = list(NCBIXML.parse(result_handle))

#writes the top 10 hits to a log file
with open("blast_top_hits.log", "w") as log_file:
    # Write the header row
    log_file.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

    # loop through each records in the BLAST results
    for record in blast_records:
        for alignment in record.alignments[:10]:  # Limit to top 10 hits and loops through the top 10 alignments for each record
            hsp = alignment.hsps[0]  # Only keep the best high scoring segment pair 

            # Extract the information and format the output
            output_line = f"{alignment.accession}\t{hsp.identities * 100 / hsp.align_length:.2f}\t{hsp.align_length}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{hsp.bits}\t{hsp.expect}\t{alignment.title}\n"

            # Write the output line to the log file
            log_file.write(output_line)


#Output file in GitHub as "Top_10_blast_top_hits.log"
```
