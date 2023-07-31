import argparse
import concurrent.futures
import csv
import os
import subprocess
import xml.etree.ElementTree as ET
from Bio import Entrez
import shutil
import pandas as pd


def download_sra(sra_id):
    print("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    subprocess.call(prefetch, shell=True)

def extract_sra(sra_id, output_dir):
    print("Generating fastq for: " + sra_id)
    fastq_dump = f"parallel-fastq-dump --sra-id {sra_id} --threads 20 --outdir {output_dir} --skip-technical --split-files"
    subprocess.call(fastq_dump, shell=True)

    # Compress the generated FASTQ files
    print("Compressing fastq for: " + sra_id)
    gzip_cmd = f"gzip {output_dir}/{sra_id}*.fastq"
    subprocess.call(gzip_cmd, shell=True)

    # Delete the directory containing the .sra file to save disk space
    print("Deleting directory for: " + sra_id)
    shutil.rmtree(f"./{sra_id}")

def fetch_sra_metadata(sra_id):
    # Fetch metadata for each file using esearch and efetch
    metadata_fetch = f'esearch -db sra -query "{sra_id}" | efetch -format runinfo'
    print("Fetching SRA metadata for: " + sra_id)
    result = subprocess.check_output(metadata_fetch, shell=True)

    # Parse the metadata and return it as a DataFrame
    data = [row.split(',') for row in result.decode('utf-8').split('\n')]
    df = pd.DataFrame(data[1:], columns=data[0])
    df.set_index('Run', inplace=True)
    return df

def fetch_sample_metadata(sra_id):
    # Set your email address (required by NCBI)
    Entrez.email = 'tomdemeter22@gmail.com'

    # Retrieve metadata for SRR accession
    handle = Entrez.esearch(db='sra', term=sra_id)
    record = Entrez.read(handle)
    handle = Entrez.efetch(db='sra', id=record['IdList'][0], rettype='full', retmode='xml')
    xml_data = handle.read()

    # Parse the XML data using ElementTree
    root = ET.fromstring(xml_data)

    # Extract information from the metadata and write it to a DataFrame
    data = {'SRR': sra_id}
    for attribute in root.findall('.//SAMPLE_ATTRIBUTE'):
        tag = attribute.find('TAG').text
        value = attribute.find('VALUE').text
        data[tag] = value
    df = pd.DataFrame(data, index=[0])
    return df

parser = argparse.ArgumentParser(description='Process SRA files.')
parser.add_argument('--input', type=str, help='Location of accession_numbers.txt')
parser.add_argument('--output', type=str, help='Output directory')
args = parser.parse_args()

# Read SRA numbers from input file
with open(args.input, 'r') as f:
    sra_numbers = f.read().splitlines()

# Create output directory if it doesn't exist
if not os.path.exists(args.output):
    os.makedirs(args.output)

# Download the .sra files in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(download_sra, sra_numbers)

# Extract the .sra files into the specified output directory one at a time
for sra_id in sra_numbers:
    extract_sra(sra_id, args.output)

# Create an empty DataFrame to store the final metadata
srr_metadata = pd.DataFrame()

# Fetch SRA and sample metadata one sample at a time and append it to the final DataFrame
for sra_id in sra_numbers:
    sample_df = fetch_sample_metadata(sra_id)
    srr_df = fetch_sra_metadata(sra_id)
  
    # Join the two dataframes on the first column of srr_df and column SRR of sample_df
    merged_df = pd.merge(sample_df, srr_df, right_index=True, left_on='SRR')
    
    # Append the merged dataframe to the final dataframe
    srr_metadata = pd.concat([srr_metadata, merged_df])

# Save the final dataframe to a CSV file named SRR_metadata.csv
srr_metadata.to_csv(f'{args.output}/SRR_metadata.csv', index=False)
