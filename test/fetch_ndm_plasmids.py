#!/usr/bin/env python3
"""
Pipeline to reproduce Figure 7 from the NDM surveillance paper.

Reference:  NDM151 (E. coli) IncHI2/IncHI2A plasmid  [SAMEA118261226]
Queries (inner → outer rings, chronological order):
    NDM152 (C. freundii)       SAMEA118261227  Jan 2023
    NDM109 (En. hormaechei)    SAMEA118261186  Feb 2023
    NDM153 (E. coli)           SAMEA118261228  Mar 2023
    NDM156 (E. coli)           SAMEA118261230  Apr 2023
    NDM144 (K. pneumoniae)     SAMEA118261219  May 2023
    NDM146 (K. pneumoniae)     SAMEA118261221  May 2023
    NDM101 (En. hormaechei)    SAMEA118261178  Jun 2023
    NDM116 (En. hormaechei)    SAMEA118261192  Jun 2023
    NDM139 (K. pneumoniae)     SAMEA118261214  Jul 2023
    NDM159 (C. brakii)         SAMEA118261233  Aug 2023
    NDM117 (En. hormaechei)    SAMEA118261193  Oct 2023
    NDM104 (En. hormaechei)    SAMEA118261185  Dec 2023
    NDM160 (K. pneumoniae)     SAMEA118261234  Nov 2023
    NDM121 (En. hormaechei)    SAMEA118261197  Nov 2023

Steps:
  1. Download assemblies from ENA/NCBI via BioSample accession
  2. Filter contigs to IncHI2/IncHI2A plasmids (blastn vs PlasmidFinder or
     length/replicon heuristic)
  3. Write per-isolate FASTA files ready for pygenomecomp-wasm upload

Requirements:
    pip install biopython requests tqdm
    conda/apt install ncbi-datasets-cli  (or use wget fallback)
    blastn in PATH (optional, for replicon typing step)

Usage:
    python fetch_ndm_plasmids.py [--outdir plasmid_fastas] [--email your@email]
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

try:
    import requests
    from Bio import Entrez, SeqIO
    from tqdm import tqdm
except ImportError:
    sys.exit(
        "Missing dependencies. Run:\n"
        "  pip install biopython requests tqdm"
    )

# ---------------------------------------------------------------------------
# Isolate manifest
# (name, biosample, species, month_label, year)
# Order matches Figure 7 inner→outer rings
# ---------------------------------------------------------------------------
REFERENCE = {
    "name": "NDM151",
    "biosample": "SAMEA118261226",
    "species": "E. coli",
    "replicon_type": "IncHI2/IncHI2A",
    "approx_size_bp": 350_000,   # ~350 kbp visible in figure
}

QUERIES = [
    # (name,    biosample,          species,           month)
    ("NDM152", "SAMEA118261227",  "C. freundii",       "Jan 2023"),
    ("NDM109", "SAMEA118261186",  "En. hormaechei",    "Feb 2023"),
    ("NDM153", "SAMEA118261228",  "E. coli",           "Mar 2023"),
    ("NDM156", "SAMEA118261230",  "E. coli",           "Apr 2023"),
    ("NDM144", "SAMEA118261219",  "K. pneumoniae",     "May 2023"),
    ("NDM146", "SAMEA118261221",  "K. pneumoniae",     "May 2023"),
    ("NDM101", "SAMEA118261178",  "En. hormaechei",    "Jun 2023"),
    ("NDM116", "SAMEA118261192",  "En. hormaechei",    "Jun 2023"),
    ("NDM139", "SAMEA118261214",  "K. pneumoniae",     "Jul 2023"),
    ("NDM159", "SAMEA118261233",  "C. brakii",         "Aug 2023"),
    ("NDM117", "SAMEA118261193",  "En. hormaechei",    "Oct 2023"),
    ("NDM104", "SAMEA118261185",  "En. hormaechei",    "Dec 2023"),
    ("NDM160", "SAMEA118261234",  "K. pneumoniae",     "Nov 2023"),
    ("NDM121", "SAMEA118261197",  "En. hormaechei",    "Nov 2023"),
]

ALL_ISOLATES = [
    (REFERENCE["name"], REFERENCE["biosample"]),
] + [(q[0], q[1]) for q in QUERIES]

# ---------------------------------------------------------------------------
# Helper: resolve BioSample → assembly accession via ENA API
# ---------------------------------------------------------------------------

def biosample_to_assembly_ena(biosample: str) -> list[str]:
    """
    Query ENA portal API to resolve a BioSample to assembly accession(s).
    Returns list of GCA/ERZ accessions.
    """
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search"
        f"?result=assembly&query=sample_accession%3D{biosample}"
        "&fields=accession,assembly_level&format=json&limit=10"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        data = r.json()
        return [d["accession"] for d in data if d.get("accession")]
    except Exception as e:
        print(f"  [WARN] ENA lookup failed for {biosample}: {e}")
        return []


def biosample_to_assembly_ncbi(biosample: str, email: str) -> list[str]:
    """
    Fallback: use NCBI Entrez to resolve BioSample → assembly accessions.
    """
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="assembly", term=f"{biosample}[BioSample]")
        record = Entrez.read(handle)
        handle.close()
        ids = record.get("IdList", [])
        accessions = []
        for uid in ids:
            h2 = Entrez.esummary(db="assembly", id=uid)
            summary = Entrez.read(h2, validate=False)
            h2.close()
            doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
            acc = doc.get("AssemblyAccession", "")
            if acc:
                accessions.append(acc)
        return accessions
    except Exception as e:
        print(f"  [WARN] NCBI Entrez lookup failed for {biosample}: {e}")
        return []


# ---------------------------------------------------------------------------
# Helper: download assembly FASTA from ENA
# ---------------------------------------------------------------------------

def download_assembly_fasta_ena(accession: str, outpath: Path) -> bool:
    """
    Download the assembled genome FASTA from ENA FTP.
    Works for ERZ* and GCA* accessions stored in ENA.
    """
    # Try ENA FTP structure

    # ERZ style
    url = (
        f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
        "?download=true&gzip=true"
    )
    return _stream_download(url, outpath, gzipped=True)


def download_assembly_fasta_ncbi(accession: str, outpath: Path) -> bool:
    """
    Download genome FASTA from NCBI datasets CLI or direct FTP.
    """
    # Try ncbi datasets CLI first (fast, handles all GCA/GCF)
    try:
        result = subprocess.run(
            [
                "datasets", "download", "genome",
                "accession", accession,
                "--include", "genome",
                "--filename", str(outpath.with_suffix(".zip")),
                "--no-progressbar",
            ],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode == 0:
            # Unzip and extract fasta
            zip_path = outpath.with_suffix(".zip")
            subprocess.run(
                ["unzip", "-o", str(zip_path),
                 "ncbi_dataset/data/*/GCA_*.fna",
                 "ncbi_dataset/data/*/GCF_*.fna",
                 "-d", str(outpath.parent)],
                capture_output=True
            )
            # Find extracted fna and move to outpath
            dataset_dir = outpath.parent / "ncbi_dataset"
            for fna in dataset_dir.rglob("*.fna"):
                fna.rename(outpath)
                zip_path.unlink(missing_ok=True)
                import shutil
                shutil.rmtree(dataset_dir, ignore_errors=True)
                return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # Fallback: NCBI FTP via Entrez summary
    return _ncbi_ftp_download(accession, outpath)


def _ncbi_ftp_download(accession: str, outpath: Path) -> bool:
    """Download genome FASTA via NCBI FTP URL from assembly summary."""
    try:
        # Build FTP path from GCA accession
        # GCA_XXXXXXXXX.N → /genomes/all/GCA/XXX/XXX/XXX/GCA_XXXXXXXXX.N_*/
        m = re.match(r"(GCA|GCF)_(\d{3})(\d{3})(\d{3})\.(\d+)", accession)
        if not m:
            return False
        prefix, p1, p2, p3, ver = m.groups()
        ftp_base = (
            f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/"
            f"{p1}/{p2}/{p3}/"
        )
        # List directory
        r = requests.get(ftp_base, timeout=20)
        # Find assembly directory
        dirs = re.findall(rf'href="({accession}[^"]+)/"', r.text)
        if not dirs:
            return False
        asm_dir = dirs[0]
        fna_url = (
            f"{ftp_base}{asm_dir}/{asm_dir}_genomic.fna.gz"
        )
        return _stream_download(fna_url, outpath, gzipped=True)
    except Exception as e:
        print(f"  [WARN] FTP download failed for {accession}: {e}")
        return False


def _stream_download(url: str, outpath: Path, gzipped: bool = False) -> bool:
    try:
        r = requests.get(url, stream=True, timeout=120)
        r.raise_for_status()
        if gzipped:
            # Download to a temporary file with .gz extension for gunzip
            temp_gzipped_path = outpath.parent / (outpath.name + ".gz.tmp")
            with open(temp_gzipped_path, "wb") as f:
                for chunk in r.iter_content(65536):
                    f.write(chunk)
            # Decompress the gzipped file content directly to the final outpath
            with open(outpath, "wb") as f_out:
                subprocess.run(["gunzip", "-c", str(temp_gzipped_path)], stdout=f_out, check=True)
            temp_gzipped_path.unlink() # Delete the temporary gzipped file
        else:
            # Download directly to a temporary file and then rename
            temp_path = outpath.with_suffix(".tmp")
            with open(temp_path, "wb") as f:
                for chunk in r.iter_content(65536):
                    f.write(chunk)
            temp_path.rename(outpath)
        return True
    except Exception as e:
        print(f"  [WARN] Download failed from {url}: {e}")
        return False


# ---------------------------------------------------------------------------
# Helper: extract IncHI2 plasmid contig(s) from assembly
# Strategy 1: size filter (IncHI2 plasmids ~200–400 kb)
# Strategy 2: blastn vs blaNDM-1 + IncHI2 rep gene
# ---------------------------------------------------------------------------

INCHI2_SIZE_MIN = 150_000   # bp
INCHI2_SIZE_MAX = 500_000   # bp

# Minimal blaNDM-1 sequence for screening (first 200 bp is sufficient as seed)
BLANDM1_SEED = """>blaNDM-1_seed
ATGAAATTCTTCGGCACCTTCGGAATGGCCGATCAGCTTGTTGAACAGGGCCTGCAAGA
GCAGGTCGCGCTGCTGGGTGACGGCTTCCAGCAGGTCTGGCACGACCAGCCCATGATCG
ATGGCTTCCAGCAAATGCAGCAAACCCAGCAGGCGGTGCTGGTGCAGGATGCGTTTGCC
GCGCTGCGTGAACGCGGTCTGCATCTGGAGCAGGCGCGCGAACAGGCGCTGCTGGGTG
"""

INCHI2_REP_SEED = """>IncHI2_rep_seed
ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCC
TGTTTTTGCTCACCCAGAAACGCTGGTAAAGTAAAGATGCTGAAGATCAGTTGGGTGCA
"""


def find_inchi2_plasmid_contigs(
    fasta_path: Path,
    isolate_name: str,
    use_blast: bool = True,
) -> list:
    """
    Return SeqRecord(s) that are likely the IncHI2/IncHI2A plasmid.
    Priority:
      1. blastn hit against blaNDM-1 AND IncHI2 rep seed on same contig
      2. blastn hit against blaNDM-1 on a contig in size range
      3. Size-only filter (150–500 kb contigs)
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print(f"  [WARN] {isolate_name}: empty FASTA")
        return []

    # ---- Strategy 1/2: blastn ----
    if use_blast and _blast_available():
        hits_ndm = _blast_seed(BLANDM1_SEED, fasta_path, "blaNDM1_seed.fa")
        hits_rep = _blast_seed(INCHI2_REP_SEED, fasta_path, "IncHI2rep_seed.fa")

        # Contigs with NDM hit + rep hit → strongest evidence
        both = set(hits_ndm) & set(hits_rep)
        if both:
            selected = [r for r in records if r.id in both]
            print(f"  {isolate_name}: {len(selected)} contig(s) — NDM+IncHI2 rep hit")
            return selected

        # NDM hit alone
        if hits_ndm:
            selected = [r for r in records if r.id in hits_ndm]
            print(f"  {isolate_name}: {len(selected)} contig(s) — blaNDM-1 hit only")
            return selected

    # ---- Strategy 3: size filter ----
    sized = [
        r for r in records
        if INCHI2_SIZE_MIN <= len(r.seq) <= INCHI2_SIZE_MAX
    ]
    if sized:
        print(
            f"  {isolate_name}: {len(sized)} contig(s) — "
            f"size filter ({INCHI2_SIZE_MIN//1000}–{INCHI2_SIZE_MAX//1000} kb)"
        )
        return sized

    # Last resort: largest contig
    biggest = sorted(records, key=lambda r: len(r.seq), reverse=True)[0]
    print(
        f"  [WARN] {isolate_name}: no IncHI2 evidence found — "
        f"using largest contig ({len(biggest.seq)//1000} kb)"
    )
    return [biggest]


def _blast_available() -> bool:
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def _blast_seed(seed_seq: str, fasta_path: Path, seed_filename: str) -> set[str]:
    """Write seed, run blastn, return set of matching subject contig IDs."""
    tmp_seed = Path(f"/tmp/{seed_filename}")
    tmp_seed.write_text(seed_seq)
    tmp_out = Path(f"/tmp/blast_out_{fasta_path.stem}.txt")
    try:
        subprocess.run(
            [
                "blastn", "-query", str(tmp_seed),
                "-subject", str(fasta_path),
                "-outfmt", "6 sseqid pident length",
                "-perc_identity", "80",
                "-out", str(tmp_out),
                "-dust", "no",
            ],
            capture_output=True, check=True, timeout=120
        )
        hits = set()
        for line in tmp_out.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and int(parts[2]) >= 100:
                hits.add(parts[0])
        return hits
    except Exception:
        return set()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run(outdir: Path, email: str, skip_existing: bool = True):
    outdir.mkdir(parents=True, exist_ok=True)
    asm_dir = outdir / "assemblies"
    asm_dir.mkdir(exist_ok=True)
    plasmid_dir = outdir / "plasmid_fastas"
    plasmid_dir.mkdir(exist_ok=True)

    manifest = {}  # name → plasmid fasta path

    print("\n=== Step 1: Resolve BioSample → Assembly accession ===\n")
    accession_map = {}  # isolate_name → [accessions]
    for name, biosample in ALL_ISOLATES:
        print(f"  {name} ({biosample}) ...", end=" ", flush=True)
        accs = biosample_to_assembly_ena(biosample)
        if not accs:
            accs = biosample_to_assembly_ncbi(biosample, email)
        if accs:
            accession_map[name] = accs
            print(f"→ {accs}")
        else:
            print("→ NOT FOUND (manual download needed)")
        time.sleep(0.3)  # be polite to APIs

    # Save manifest
    (outdir / "accession_map.json").write_text(
        json.dumps(accession_map, indent=2)
    )

    print("\n=== Step 2: Download assemblies ===\n")
    for name, biosample in ALL_ISOLATES:
        accs = accession_map.get(name, [])
        if not accs:
            print(f"  [SKIP] {name}: no accession found")
            continue

        acc = accs[0]  # use first (usually best-quality) assembly
        fasta_out = asm_dir / f"{name}_{acc}_genomic.fna"

        if skip_existing and fasta_out.exists():
            print(f"  [SKIP] {name}: {fasta_out.name} already exists")
            continue

        print(f"  Downloading {name} ({acc}) ...", end=" ", flush=True)
        success = download_assembly_fasta_ena(acc, fasta_out)
        if not success:
            success = download_assembly_fasta_ncbi(acc, fasta_out)
        if success:
            print(f"OK -> {fasta_out.name}")
        else:
            print(f"FAILED - try manual download from:")
            print(f"    https://www.ebi.ac.uk/ena/browser/view/{biosample}")
        time.sleep(1)

    print("\n=== Step 3: Extract IncHI2/IncHI2A plasmid contig(s) ===\n")
    for name, biosample in ALL_ISOLATES:
        accs = accession_map.get(name, [])
        if not accs:
            continue
        acc = accs[0]
        fasta_in = asm_dir / f"{name}_{acc}_genomic.fna"
        if not fasta_in.exists():
            print(f"  [SKIP] {name}: assembly not downloaded")
            continue

        plasmid_fasta = plasmid_dir / f"{name}_IncHI2_plasmid.fna"
        if skip_existing and plasmid_fasta.exists():
            print(f"  [SKIP] {name}: plasmid FASTA already exists")
            manifest[name] = plasmid_fasta
            continue

        contigs = find_inchi2_plasmid_contigs(fasta_in, name, use_blast=True)
        if contigs:
            # Rename record for clarity
            for i, rec in enumerate(contigs):
                rec.id = f"{name}_IncHI2_contig{i+1}"
                rec.description = f"{name} IncHI2/IncHI2A plasmid contig {i+1}"
            SeqIO.write(contigs, plasmid_fasta, "fasta")
            manifest[name] = plasmid_fasta

    print("\n=== Step 4: Summary — files ready for pygenomecomp-wasm ===\n")
    ref_path = manifest.get(REFERENCE["name"])
    query_paths = [manifest.get(q[0]) for q in QUERIES]

    print(f"  Reference FASTA (upload as 'Reference Genome'):")
    if ref_path and ref_path.exists():
        size_kb = ref_path.stat().st_size // 1024
        print(f"    ✓ {ref_path}  ({size_kb} KB)")
    else:
        print(f"    ✗ NDM151 plasmid FASTA not generated — check download")

    print(f"\n  Query FASMTAs (upload as 'Query Genomes', in this order = inner→outer):")
    for i, (q, path) in enumerate(zip(QUERIES, query_paths)):
        name, biosample, species, month = q
        if path and path.exists():
            size_kb = path.stat().st_size // 1024
            print(f"    Ring {i+1:2d}  ✓  {path.name}  — {name} ({species}) {month}")
        else:
            print(f"    Ring {i+1:2d}  ✗  {name} ({species}) {month} — missing")

    print("\n  Recommended pygenomecomp-wasm settings:")
    print("    Mode:              Circular")
    print("    Min Identity:      70%")
    print("    Min Coverage:      0%")
    print("    Min Alignment Len: 500 bp")
    print("\n  Then upload the reference + all query FASTAs, click 'Run Comparison',")
    print("  and download the SVG.\n")

    # Write a combined multi-FASTA of all queries (optional convenience)
    combined = plasmid_dir / "ALL_QUERIES_IncHI2_plasmids.fna"
    all_q_records = []
    for q, path in zip(QUERIES, query_paths):
        if path and path.exists():
            all_q_records.extend(SeqIO.parse(path, "fasta"))
    if all_q_records:
        SeqIO.write(all_q_records, combined, "fasta")
        print(f"  Combined query FASTA (all queries in one file):")
        print(f"    {combined}\n")

    return manifest


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch NDM IncHI2 plasmid sequences for pygenomecomp-wasm Figure 7"
    )
    parser.add_argument(
        "--outdir", default="ndm_plasmid_figure7",
        help="Output directory (default: ndm_plasmid_figure7)"
    )
    parser.add_argument(
        "--email", default="motroy@post.bgu.ac.il",
        help="Email for NCBI Entrez (required by NCBI policy)"
    )
    parser.add_argument(
        "--no-skip", action="store_true",
        help="Re-download even if files already exist"
    )
    args = parser.parse_args()

    run(
        outdir=Path(args.outdir),
        email=args.email,
        skip_existing=not args.no_skip,
    )
def find_inchi2_plasmid_contigs(
    fasta_path: Path,
    isolate_name: str,
    use_blast: bool = True,
) -> list:
    """
    Return SeqRecord(s) that are likely the IncHI2/IncHI2A plasmid.
    Priority:
      1. blastn hit against blaNDM-1 AND IncHI2 rep seed on same contig
      2. blastn hit against blaNDM-1 on a contig in size range
      3. Size-only filter (150–500 kb contigs)
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        print(f"  [WARN] {isolate_name}: empty FASTA")
        return []

    # ---- Strategy 1/2: blastn ----
    if use_blast and _blast_available():
        hits_ndm = _blast_seed(BLANDM1_SEED, fasta_path, "blaNDM1_seed.fa")
        hits_rep = _blast_seed(INCHI2_REP_SEED, fasta_path, "IncHI2rep_seed.fa")

        # Contigs with NDM hit + rep hit → strongest evidence
        both = set(hits_ndm) & set(hits_rep)
        if both:
            selected = [r for r in records if r.id in both]
            print(f"  {isolate_name}: {len(selected)} contig(s) — NDM+IncHI2 rep hit")
            return selected

        # NDM hit alone
        if hits_ndm:
            selected = [r for r in records if r.id in hits_ndm]
            print(f"  {isolate_name}: {len(selected)} contig(s) — blaNDM-1 hit only")
            return selected

    # ---- Strategy 3: size filter ----
    sized = [
        r for r in records
        if INCHI2_SIZE_MIN <= len(r.seq) <= INCHI2_SIZE_MAX
    ]
    if sized:
        print(
            f"  {isolate_name}: {len(sized)} contig(s) — "
            f"size filter ({INCHI2_SIZE_MIN//1000}–{INCHI2_SIZE_MAX//1000} kb)"
        )
        return sized

    # Last resort: largest contig
    biggest = sorted(records, key=lambda r: len(r.seq), reverse=True)[0]
    print(
        f"  [WARN] {isolate_name}: no IncHI2 evidence found — "
        f"using largest contig ({len(biggest.seq)//1000} kb)"
    )
    return [biggest]


def _blast_available() -> bool:
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def _blast_seed(seed_seq: str, fasta_path: Path, seed_filename: str) -> set[str]:
    """Write seed, run blastn, return set of matching subject contig IDs."""
    tmp_seed = Path(f"/tmp/{seed_filename}")
    tmp_seed.write_text(seed_seq)
    tmp_out = Path(f"/tmp/blast_out_{fasta_path.stem}.txt")
    try:
        subprocess.run(
            [
                "blastn", "-query", str(tmp_seed),
                "-subject", str(fasta_path),
                "-outfmt", "6 sseqid pident length",
                "-perc_identity", "80",
                "-out", str(tmp_out),
                "-dust", "no",
            ],
            capture_output=True, check=True, timeout=120
        )
        hits = set()
        for line in tmp_out.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and int(parts[2]) >= 100:
                hits.add(parts[0])
        return hits
    except Exception:
        return set()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run(outdir: Path, email: str, skip_existing: bool = True):
    outdir.mkdir(parents=True, exist_ok=True)
    asm_dir = outdir / "assemblies"
    asm_dir.mkdir(exist_ok=True)
    plasmid_dir = outdir / "plasmid_fastas"
    plasmid_dir.mkdir(exist_ok=True)

    manifest = {}  # name → plasmid fasta path

    print("\n=== Step 1: Resolve BioSample → Assembly accession ===\n")
    accession_map = {}  # isolate_name → [accessions]
    for name, biosample in ALL_ISOLATES:
        print(f"  {name} ({biosample}) ...", end=" ", flush=True)
        accs = biosample_to_assembly_ena(biosample)
        if not accs:
            accs = biosample_to_assembly_ncbi(biosample, email)
        if accs:
            accession_map[name] = accs
            print(f"→ {accs}")
        else:
            print("→ NOT FOUND (manual download needed)")
        time.sleep(0.3)  # be polite to APIs

    # Save manifest
    (outdir / "accession_map.json").write_text(
        json.dumps(accession_map, indent=2)
    )

    print("\n=== Step 2: Download assemblies ===\n")
    for name, biosample in ALL_ISOLATES:
        accs = accession_map.get(name, [])
        if not accs:
            print(f"  [SKIP] {name}: no accession found")
            continue

        acc = accs[0]  # use first (usually best-quality) assembly
        fasta_out = asm_dir / f"{name}_{acc}_genomic.fna"

        if skip_existing and fasta_out.exists():
            print(f"  [SKIP] {name}: {fasta_out.name} already exists")
            continue

        print(f"  Downloading {name} ({acc}) ...", end=" ", flush=True)
        success = download_assembly_fasta_ena(acc, fasta_out)
        if not success:
            success = download_assembly_fasta_ncbi(acc, fasta_out)
        if success:
            print(f"OK -> {fasta_out.name}")
        else:
            print(f"FAILED - try manual download from:")
            print(f"    https://www.ebi.ac.uk/ena/browser/view/{biosample}")
        time.sleep(1)

    print("\n=== Step 3: Extract IncHI2/IncHI2A plasmid contig(s) ===\n")
    for name, biosample in ALL_ISOLATES:
        accs = accession_map.get(name, [])
        if not accs:
            continue
        acc = accs[0]
        fasta_in = asm_dir / f"{name}_{acc}_genomic.fna"
        if not fasta_in.exists():
            print(f"  [SKIP] {name}: assembly not downloaded")
            continue

        plasmid_fasta = plasmid_dir / f"{name}_IncHI2_plasmid.fna"
        if skip_existing and plasmid_fasta.exists():
            print(f"  [SKIP] {name}: plasmid FASTA already exists")
            manifest[name] = plasmid_fasta
            continue

        contigs = find_inchi2_plasmid_contigs(fasta_in, name, use_blast=True)
        if contigs:
            # Rename record for clarity
            for i, rec in enumerate(contigs):
                rec.id = f"{name}_IncHI2_contig{i+1}"
                rec.description = f"{name} IncHI2/IncHI2A plasmid contig {i+1}"
            SeqIO.write(contigs, plasmid_fasta, "fasta")
            manifest[name] = plasmid_fasta

    print("\n=== Step 4: Summary — files ready for pygenomecomp-wasm ===\n")
    ref_path = manifest.get(REFERENCE["name"])
    query_paths = [manifest.get(q[0]) for q in QUERIES]

    print(f"  Reference FASTA (upload as 'Reference Genome'):")
    if ref_path and ref_path.exists():
        size_kb = ref_path.stat().st_size // 1024
        print(f"    ✓ {ref_path}  ({size_kb} KB)")
    else:
        print(f"    ✗ NDM151 plasmid FASTA not generated — check download")

    print(f"\n  Query FASMTAs (upload as 'Query Genomes', in this order = inner→outer):")
    for i, (q, path) in enumerate(zip(QUERIES, query_paths)):
        name, biosample, species, month = q
        if path and path.exists():
            size_kb = path.stat().st_size // 1024
            print(f"    Ring {i+1:2d}  ✓  {path.name}  — {name} ({species}) {month}")
        else:
            print(f"    Ring {i+1:2d}  ✗  {name} ({species}) {month} — missing")

    print("\n  Recommended pygenomecomp-wasm settings:")
    print("    Mode:              Circular")
    print("    Min Identity:      70%")
    print("    Min Coverage:      0%")
    print("    Min Alignment Len: 500 bp")
    print("\n  Then upload the reference + all query FASTAs, click 'Run Comparison',")
    print("  and download the SVG.\n")

    # Write a combined multi-FASTA of all queries (optional convenience)
    combined = plasmid_dir / "ALL_QUERIES_IncHI2_plasmids.fna"
    all_q_records = []
    for q, path in zip(QUERIES, query_paths):
        if path and path.exists():
            all_q_records.extend(SeqIO.parse(path, "fasta"))
    if all_q_records:
        SeqIO.write(all_q_records, combined, "fasta")
        print(f"  Combined query FASTA (all queries in one file):")
        print(f"    {combined}\n")

    return manifest


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch NDM IncHI2 plasmid sequences for pygenomecomp-wasm Figure 7"
    )
    parser.add_argument(
        "--outdir", default="ndm_plasmid_figure7",
        help="Output directory (default: ndm_plasmid_figure7)"
    )
    parser.add_argument(
        "--email", default="motroy@post.bgu.ac.il",
        help="Email for NCBI Entrez (required by NCBI policy)"
    )
    parser.add_argument(
        "--no-skip", action="store_true",
        help="Re-download even if files already exist"
    )
    args = parser.parse_args()

    run(
        outdir=Path(args.outdir),
        email=args.email,
        skip_existing=not args.no_skip,
    )
