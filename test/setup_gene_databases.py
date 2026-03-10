import os
from Bio import SeqIO
import requests

def download_file(url, filename):
    if not os.path.exists(filename):
        print(f'[INFO] Downloading {filename}...')
        response = requests.get(url)
        with open(filename, 'wb') as f:
            f.write(response.content)
    else:
        print(f'[SKIP] {filename} already exists.')

def create_subsets():
    # URLs for Abricate databases
    plasmidfinder_url = 'https://raw.githubusercontent.com/tseemann/abricate/refs/heads/master/db/plasmidfinder/sequences'
    ncbi_url = 'https://raw.githubusercontent.com/tseemann/abricate/refs/heads/master/db/ncbi/sequences'

    # Step 1: Download raw databases
    download_file(plasmidfinder_url, 'plasmidfinder.seqs.fasta')
    download_file(ncbi_url, 'ncbi.seqs.fasta')

    # Step 2: Subset NDM variants
    print('[INFO] Extracting blaNDM variants...')
    ndm_variants = []
    for record in SeqIO.parse('ncbi.seqs.fasta', 'fasta'):
        if 'blaNDM' in record.id or 'blaNDM' in record.description:
            ndm_variants.append(record)
    
    with open('blaNDM_variants.fasta', 'w') as f:
        SeqIO.write(ndm_variants, f, 'fasta')
    print(f'[SUCCESS] Created blaNDM_variants.fasta ({len(ndm_variants)} sequences)')

    # Step 3: Subset IncHI replicons
    print('[INFO] Extracting IncHI replicons...')
    inchi_replicons = []
    for record in SeqIO.parse('plasmidfinder.seqs.fasta', 'fasta'):
        if 'IncHI' in record.id or 'IncHI' in record.description:
            inchi_replicons.append(record)
    
    with open('IncHI_replicons.fasta', 'w') as f:
        SeqIO.write(inchi_replicons, f, 'fasta')
    print(f'[SUCCESS] Created IncHI_replicons.fasta ({len(inchi_replicons)} sequences)')

if __name__ == '__main__':
    create_subsets()
