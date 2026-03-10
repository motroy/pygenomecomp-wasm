#!/usr/bin/env python3
import os
import subprocess
import shutil
import requests
from pathlib import Path
from Bio import SeqIO

BLANDM_DB = 'blaNDM_variants.fasta'
INCHI_DB = 'IncHI_replicons.fasta'

REFERENCE = {'name': 'NDM151', 'biosample': 'SAMEA118261226', 'species': 'E. coli'}
QUERIES = [
    ('NDM152', 'SAMEA118261227', 'C. freundii', 'Jan 2023'),
    ('NDM109', 'SAMEA118261186', 'En. hormaechei', 'Feb 2023'),
    ('NDM153', 'SAMEA118261228', 'E. coli', 'Mar 2023'),
    ('NDM156', 'SAMEA118261230', 'E. coli', 'Apr 2023'),
    ('NDM144', 'SAMEA118261219', 'K. pneumoniae', 'May 2023'),
    ('NDM146', 'SAMEA118261221', 'K. pneumoniae', 'May 2023'),
    ('NDM101', 'SAMEA118261178', 'En. hormaechei', 'Jun 2023'),
    ('NDM116', 'SAMEA118261192', 'En. hormaechei', 'Jun 2023'),
    ('NDM139', 'SAMEA118261214', 'K. pneumoniae', 'Jul 2023'),
    ('NDM159', 'SAMEA118261233', 'C. brakii', 'Aug 2023'),
    ('NDM117', 'SAMEA118261193', 'En. hormaechei', 'Oct 2023'),
    ('NDM104', 'SAMEA118261185', 'En. hormaechei', 'Dec 2023'),
    ('NDM160', 'SAMEA118261234', 'K. pneumoniae', 'Nov 2023'),
    ('NDM121', 'SAMEA118261197', 'En. hormaechei', 'Nov 2023'),
]

def _blast_db(db_path, query_fasta):
    out_fmt = '6 qseqid sseqid bitscore'
    cmd = ['blastn', '-query', query_fasta, '-subject', db_path, '-outfmt', out_fmt]
    hits = {}
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in res.stdout.strip().split('\n'):
            if not line: continue
            qsid, _, bitscore = line.split('\t')
            hits[qsid] = max(hits.get(qsid, 0), float(bitscore))
    except: pass
    return hits

def download_assembly(biosample, out_path):
    # Using Entrez to find accession and NCBI Datasets for the download
    try:
        # First, search via Entrez to get assembly accession
        import requests
        search_url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/biosample/{biosample}/dataset_report'
        resp = requests.get(search_url).json()
        accession = resp['reports'][0]['accession']
        
        # Use wget to download the assembly package
        zip_file = f"{biosample}.zip"
        dl_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA"
        subprocess.run(["wget", "-q", dl_url, "-O", zip_file], check=True)
        
        # Unzip and find the .fna file
        subprocess.run(["unzip", "-o", zip_file, "-d", biosample], capture_output=True)
        for fna in Path(biosample).rglob("*.fna"):
            shutil.copy(str(fna), out_path)
            break
        # Clean up
        shutil.rmtree(biosample)
        os.remove(zip_file)
        return os.path.exists(out_path)
    except:
        return False


def process_isolate(name, biosample, species, outdir, all_fna):
    print(f'  Processing {name} ({biosample})...')
    tmp_fna = os.path.join(outdir, f"{name}_temp.fna")
    if not download_assembly(biosample, tmp_fna):
        print(f'    [ERROR] Download failed for {name}')
        return

    ndm_hits = _blast_db(BLANDM_DB, tmp_fna)
    inchi_hits = _blast_db(INCHI_DB, tmp_fna)

    best_c = None
    max_s = -1
    all_records = list(SeqIO.parse(tmp_fna, 'fasta'))
    
    for rec in all_records:
        if rec.id in ndm_hits and rec.id in inchi_hits:
            s = ndm_hits[rec.id] + inchi_hits[rec.id]
            if 150000 <= len(rec.seq) <= 500000: s += 1000
            if s > max_s: max_s, best_c = s, rec

    if best_c:
        best_c.id = f"{name}_IncHI2_plasmid"
        SeqIO.write(best_c, os.path.join(outdir, f"{name}_IncHI2_plasmid.fna"), "fasta")
        with open(all_fna, 'a') as f: SeqIO.write(best_c, f, 'fasta')
        print(f'    [SUCCESS] Isolated contig for {name}')
    else:
        # Diagnostic Logging Section
        ndm_contigs = [rec for rec in all_records if rec.id in ndm_hits]
        inchi_contigs = [rec for rec in all_records if rec.id in inchi_hits]
        
        print(f'    [DIAGNOSTIC] {name} summary:')
        print(f'      - Total unique contigs with blaNDM: {len(ndm_contigs)}')
        print(f'      - Total unique contigs with IncHI: {len(inchi_contigs)}')
        
        all_hit_ids = set(ndm_hits.keys()) | set(inchi_hits.keys())
        if all_hit_ids:
            print(f'      - Contig lengths with hits:')
            for rec in all_records:
                if rec.id in all_hit_ids:
                    marker_type = []
                    if rec.id in ndm_hits: marker_type.append("blaNDM")
                    if rec.id in inchi_hits: marker_type.append("IncHI")
                    print(f'        * {rec.id}: {len(rec.seq)} bp (Markers: {", ".join(marker_type)})')
        else:
            print(f'      - No marker hits found on any contigs.')
        
        print(f'    [FAILED] No single contig found with both markers for {name}')

    if os.path.exists(tmp_fna): os.remove(tmp_fna)


if __name__ == "__main__":
    out = 'ndm_plasmid_figure7/plasmid_fastas'
    os.makedirs(out, exist_ok=True)
    all_fna = os.path.join(out, 'ALL_QUERIES_IncHI2_plasmids.fna')
    if os.path.exists(all_fna): os.remove(all_fna)
    
    process_isolate(REFERENCE['name'], REFERENCE['biosample'], REFERENCE['species'], out, all_fna)
    for q in QUERIES:
        process_isolate(q[0], q[1], q[2], out, all_fna)
    
    if os.path.exists(all_fna):
        subprocess.run(['gzip', '-f', '-k', all_fna])
