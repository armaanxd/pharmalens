#!/bin/bash
# PharmaLens Installation Script for Mac M3
# Run with: bash install.sh

set -e  # Exit on error

echo "🧬 PharmaLens Installation for Mac M3"
echo "======================================"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Anaconda or Miniconda first."
    exit 1
fi

# Activate bioai environment
echo -e "\n📦 Activating bioai environment..."
eval "$(conda shell.bash hook)"
conda activate bioai

# Install conda packages
echo -e "\n📦 Installing RDKit and Vina via conda..."
conda install -c conda-forge rdkit vina -y

# Install pip packages
echo -e "\n📦 Installing Python packages..."
pip install -r requirements.txt

# Download P2Rank (optional)
echo -e "\n🔍 Downloading P2Rank for binding site prediction..."
if [ ! -d "tools/p2rank_2.4.1" ]; then
    mkdir -p tools
    cd tools
    curl -L https://github.com/rdk/p2rank/releases/download/2.4.1/p2rank_2.4.1.tar.gz -o p2rank.tar.gz
    tar -xzf p2rank.tar.gz
    rm p2rank.tar.gz
    cd ..
    echo "✓ P2Rank installed to tools/p2rank_2.4.1/"
else
    echo "✓ P2Rank already installed"
fi

# Verify installation
echo -e "\n🧪 Verifying installation..."
python -c "
import rdkit
import torch
import transformers
from Bio import SeqIO
import py3Dmol
import gradio as gr

print('✅ All core packages installed successfully!')
print(f'✅ PyTorch version: {torch.__version__}')
print(f'✅ MPS (M3 GPU) available: {torch.backends.mps.is_available()}')
"

# Check Vina
if command -v vina &> /dev/null; then
    echo "✅ AutoDock Vina installed: $(vina --version 2>&1 | head -n1)"
else
    echo "⚠️  Vina not found in PATH (may need manual setup)"
fi

echo -e "\n✅ Installation complete!"
echo -e "\nNext steps:"
echo "  1. Test protein processor: python src/protein_processor.py"
echo "  2. Start building the docking engine"
echo "  3. Run the full app when ready: python app.py"
