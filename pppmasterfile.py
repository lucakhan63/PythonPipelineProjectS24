# Importing necessary modules
from Bio import Entrez, SeqIO
import os
import statistics
from pathlib import Path
import argparse

# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Track 1 Differential Expression Pipeline')
    parser.add_argument('output_dir', help='Output directory path', type=str)
    return parser.parse_args()

# Function to fetch the HCMV genome using its accession ID
def get_HCMVgenome(accessionID):
    # Setting up email for Entrez (required by NCBI)
    Entrez.email = "akhan63@luc.edu"
    #place your entrez address here
    # Fetching the HCMV genome record from NCBI
    handle = Entrez.efetch(db="nucleotide", id=accessionID, rettype="gb", retmode="text")
    # Parsing the GenBank record using BioPython's SeqIO module
    record = SeqIO.read(handle, "genbank")
    return record

# Function to extract coding sequences (CDS) from the genome record and write them to a file
def get_CDS(record, outfile):
    with open(outfile, "w") as file:
        CDS = 0
        # Iterating through features in the genome record
        for feature in record.features:
            # Checking if the feature is a CDS
            if feature.type == "CDS":
                CDS += 1
                # Extracting protein ID and sequence for the CDS
                proteinID = feature.qualifiers["protein_id"][0]
                sequence = feature.extract(record.seq)
                # Writing CDS sequence to file in FASTA format
                SeqIO.write(sequence, file, "fasta")
    return CDS

# Function to create kallisto index for transcriptome quantification
def transcriptome_index(fasta, index):
    index_str = str(index)
    fasta_str = str(fasta)
    # Running kallisto index command using os.system
    Indexcommand = "kallisto index -i " + index_str + " " + fasta_str
    os.system(Indexcommand)

# Function to quantify gene expression using kallisto
def TPM_quantification(index, sampleID, fastq1, fastq2, output_dir):
    output_directory = output_dir / f"TPM_Quantification_{sampleID}"
    output_directory.mkdir(parents=True, exist_ok=True)
    # Running kallisto quant command using os.system
    TPMcommand = f"kallisto quant -i {index} -o {output_dir} -b 30 -t 4 {fastq1} {fastq2}"
    os.system(TPMcommand)
    # Moving abundance.tsv file to sample directory
    abundance_file = output_dir / "abundance.tsv"
    new_abundance_file = output_directory / "abundance.tsv"
    os.rename(abundance_file, new_abundance_file)

# Function to calculate TPM statistics
def TPM_statistics(output_dir):
    abundance = output_dir / "abundance.tsv"
    TPM_values = []
    # Reading abundance.tsv file and extracting TPM values
    with open(abundance, "r") as file:
        next(file)
        for line in file:
            fields = line.strip().split("\t")
            TPM_values.append(float(fields[4]))
    # Calculating statistical values
    TPM_min = min(TPM_values)
    TPM_median = statistics.median(TPM_values)
    TPM_mean = statistics.mean(TPM_values)
    TPM_max = max(TPM_values)
    return str(round(TPM_min, 3)), str(round(TPM_median, 3)), str(round(TPM_mean, 3)), str(round(TPM_max, 3))

# Function to process SRR samples for quantification
def process_SRRs(samples, output_dir, index, fastqs):
    for sampleID, sample_info in samples.items():
        condition = sample_info["condition"]
        fastq1, fastq2 = sample_info["fastq_pair"]
        # Creating output directory for each sample
        TPM_directory = output_dir / f"TPM_Quantification_{sampleID}"
        TPM_directory.mkdir(parents=True, exist_ok=True)
        # Running TPM quantification
        TPM_quantification(index, sampleID, fastq1, fastq2, output_dir)
        # Calculating TPM statistics
        min_tpm, med_tpm, mean_tpm, max_tpm = TPM_statistics(TPM_directory)

# Function to get the most differentially expressed CDS
def get_mostDE_CDS(results_file):
    with open(results_file, 'r') as file:
        next(file)
        line = file.readline().strip()
        DE_CDS_id = line.split()[0]
    return DE_CDS_id

# Function to get the fasta sequence of the most differentially expressed CDS
def get_mostDE_fasta(CDS_id):
    Entrez.email = "akhan63@luc.edu"
    #place your entrez address here
    handle = Entrez.efetch(db="protein", id=CDS_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return record

# Function to download HCMV sequences from NCBI
def download_HCMVseqs(searchterm):
    search_command = "datasets download virus genome taxon " + searchterm + " --refseq --include genome"
    os.system(search_command)
    unzipcommand = "unzip ncbi_dataset.zip"
    os.system(unzipcommand)

# Function to create a local BLAST database
def create_localdb(file_input, database):
    file_input_str = str(file_input)
    database_str = str(database)
    makeblast_command = "makeblastdb -in " + file_input_str + " -out " + database_str + " -title " + database_str + " -dbtype nucl"
    os.system(makeblast_command)

# Function to run BLAST search
def run_blast(query_fasta, database, out_put):
    blast_command = "tblastn -query " + query_fasta + " -db " + database + " -out " + out_put + " -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -max_target_seqs 10"
    os.system(blast_command)

# Main function
def main():
    # Parsing command-line arguments
    args = parse_arguments()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    nuc_seqs = Path("ncbi_dataset/data/genomic.fna")
    output_file = output_dir / "PipelineProject.log"
    r_script = Path("mycode.R")
    r_results = Path("significant_transcripts.txt")

    # Fetching HCMV genome
    hcmv_acc = "NC_006273.2"
    hcmv_genome = get_HCMVgenome(hcmv_acc)

    # Getting CDS from HCMV genome
    CDS_file = output_dir / "CDS.fasta"
    CDS_count = get_CDS(hcmv_genome, CDS_file)

    # Creating kallisto index for CDS
    indexID = output_dir / "HCMV_index.idx"
    transcriptome_index(CDS_file, indexID)

    # Sample information
    samples = {
        "SRR5660030": {
            "condition": "2dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660030_1.fastq", "PipelineFastqs/SRR5660030_2.fastq")},
        "SRR5660033": {
            "condition": "6dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660033_1.fastq", "PipelineFastqs/SRR5660033_2.fastq")},
        "SRR5660044": {
            "condition": "2dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660044_1.fastq", "PipelineFastqs/SRR5660044_2.fastq")},
        "SRR5660045": {
            "condition": "6dpi",
            "fastq_pair": ("PipelineFastqs/SRR5660045_1.fastq", "PipelineFastqs/SRR5660045_2.fastq")}
    }

    # Processing SRR samples
    fastqs = sorted(Path("PipelineFastqs").glob("*.fastq"))
    process_SRRs(samples, output_dir, indexID, fastqs)

    # Writing results to log file
    with open(output_file, "a") as file:
        file.write("The HCMV genome " + hcmv_acc + " has " + str(CDS_count) + " CDS.\n")
        file.write("{:<12}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\n".format("sample", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"))
        for sampleID, condition_info in samples.items():
            condition = condition_info["condition"]
            TPM_directory = output_dir / f"TPM_Quantification_{sampleID}"
            min_tpm, med_tpm, mean_tpm, max_tpm = TPM_statistics(TPM_directory)
            file.write("{:<12}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\t{:>8}\n".format(sampleID, condition, min_tpm, med_tpm, mean_tpm, max_tpm))

    # Running R script
    os.system("Rscript " + str(r_script))

    # Writing headers for sleuth R script DE analysis values
    with open(output_file, "a") as file:
        file.write("{:<12}\t{:<20}\t{:<20}\t{:<20}\n".format("target_id", "test_stat", "pval", "qval"))

    # Writing results from R script to log file
    with open(r_results, "r") as results_file:
        next(results_file)
        for line in results_file:
            info = line.strip().split()
            target_id = info[0]
            pval = info[1]
            qval = info[2]
            test_stat = info[3]
            with open(output_file, "a") as log_file:
                log_file.write("{:<12}\t{:<20}\t{:<20}\t{:<20}\n".format(target_id, test_stat, pval, qval))

    # Getting most differentially expressed CDS
    mostDE_CDS = get_mostDE_CDS(r_results)
    DE_fasta = get_mostDE_fasta(mostDE_CDS)

    # Downloading HCMV sequences
    searchterm = "Betaherpesvirinae"
    download_HCMVseqs(searchterm)
    create_localdb(nuc_seqs, searchterm)

    # Writing most DE CDS to query fasta file
    SeqIO.write(DE_fasta, "query.fasta", "fasta")
    query = "query.fasta"
    blast_results = "HCMV_tblastn_results.csv"
    run_blast(query, searchterm, blast_results)

    # Writing top BLAST hits to log file
    with open(blast_results, 'r') as file:
        with open(output_file, 'a') as log:
            log.write("{:<12}\t{:>8}\t{:>6}\t{:>6}\t{:>6}\t{:>8}\t{:>8}\t{:>8}\t{:>10}\t{}\n".format("sacc", "pident", "length", "qstart", "qend", "sstart", "send", "bitscore", "evalue", "stitle"))
            line_count = 0
            for line in file:
                if line_count >= 10:
                    break
                values = line.strip().split('\t')
                formatted_line = "{:<12}\t{:>8}\t{:>6}\t{:>6}\t{:>6}\t{:>8}\t{:>8}\t{:>8}\t{:>10}\t{}\n".format(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9])
                log.write(formatted_line)
                line_count += 1

if __name__ == "__main__":
    main()
