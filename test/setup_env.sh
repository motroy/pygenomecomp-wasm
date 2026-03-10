#!/bin/bash
# Environment setup script for genomic analysis (Headless/Non-interactive)

# Set non-interactive frontend for Debian/Ubuntu
export DEBIAN_FRONTEND=noninteractive

echo "Installing Python dependencies..."
pip install biopython requests tqdm

echo "Updating system packages..."
apt-get update -y

echo "Installing NCBI Entrez Direct (Non-interactive)..."
# Pipe 'y' to the installer to handle the activation prompt automatically
echo "y" | sh -c "$(wget -qO- https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

echo "Updating PATH for Entrez Direct..."
export PATH=$PATH:/root/edirect

echo "Installing BLAST+..."
apt-get install -y ncbi-blast+

echo "Environment setup complete."
