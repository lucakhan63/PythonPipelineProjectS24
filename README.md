# PythonPipelineProjectS24

This is the repository for my (Ali Khan) attempt at completing the requirements of the python pipeline project for COMP 383 for the Spring 2024 semester for TRACK 1.

### Main Python Script
#### pppmasterfile.py
##### To run this script, please clone this github repository with:
```
git clone https://github.com/lucakhan63/PythonPipelineProjectS24.git
```
The SRR IDs for the test data were: SRR5660030 SRR5660033 SRR5660044 SRR5660045 

##### To run the script:
Argument format: 
```
python pppmasterfile.py <input_directory> <output_directory_file_path>
```

**Package Dependencies**

- os (used to run command line commands from Python)
 - sratoolkit 
 - kallisto 
 - statistics 
 - sleuth 
 - dplyr 
 - blast+ 
 - BioPython  
 - pathlib (used for finding, establishing, and accessing files and directories from the user's system with minimal user input)
 - argparse


**Step 1 Directions:**

- The SRAtoolkit package was installed for the purposes of converting SRA files into paired-end fastq reads
- SRA accession/SRR numbers (in the format of SRR____) for each sample were collected by visiting the respective NCBI page and looking under the "Runs" subheader, pages are listed below:
```
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
```
- The prefetch SRA-Toolkit function was utilized followed by the SRA accession/SRR number, downloading the SRA file to user's current directory. 
- fasterq-dump SRA-Toolkit function followed by the SRA accession/SRR number to convert SRR file to a paired0end fastq file. The --split-files & --skip-technicals options were used to create separate files for each paired-end read and skip technical reads respectively. After executing this function, two fastq files will be present for each sample in the user's respective, specified, directory.


**Objectives of Project**
Quantify TPM in each sample using kallisto, but first, you need to build a transcriptome index for HCMV
(NCBI accession NC_006273.2). Use Biopython to retrieve and generate the appropriate input and then build the
index with kallisto. You should extract the coding sequence (CDS) features from the GenBank format file and make a
fasta file with the RefSeq protein_id as the CDS identifier. Write the following to your log file (replace # with the
number of coding sequences in the HCMV genome):

_The HCMV genome (NC_006273.2) has # CDS._

Quantify the TPM of each CDS in each transcriptome using kallisto. For each sample (SRR number), calculate the
minimum, median, mean, and maximum TPM from the results in the abundance.tsv kallisto output file. Write
the following details to your log file, include a header row, and tab-delimit each item:

_sample condition min_tpm med_tpm mean_tpm max_tpm_

Use the output from kallisto as input for the R package sleuth to find differentially expressed genes between the
two timepoints/conditions (2dpi and 6dpi). Write the following details for each significant transcript (FDR < 0.05) to
your log file, include a header row, and tab-delimit each item:

_target_id test_stat pval qval_

What other virus strains have gene(s) encoding the most differentially expressed protein? Use Biopython to
retrieve a protein fasta file of the most differentially expressed CDS in the reference genome. Use this protein fasta
file as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily.
Think, which blast should you use? You will need to make a local database of just nucleotide sequences from the
Betaherpesvirinae subfamily. Your blast+ run should only keep the best alignment (HSP) for any single query-subject
pair of sequences. For the top 10 hits, write the following to your log file: 

_Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title._

Include the following header row in the log file, followed by the top 10 hits, and tab-delimit each item:
_sacc pident length qstart qend sstart send bitscore evalue stitle_
