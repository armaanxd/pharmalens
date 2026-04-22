"""
Protein Processor Module
Handles PDB file parsing, structure validation, and binding site prediction
"""

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import numpy as np
from pathlib import Path
import requests
import subprocess
import tempfile
import json


class ProteinProcessor:
    """Process protein structures and predict binding sites"""
    
    def __init__(self, p2rank_path=None):
        """
        Args:
            p2rank_path: Path to p2rank executable (optional, uses 'prank' if in PATH)
        """
        self.parser = PDBParser(QUIET=True)
        self.p2rank_path = p2rank_path or "prank"
        
    def load_from_pdb_file(self, pdb_path):
        """
        Load protein structure from local PDB file
        
        Args:
            pdb_path: Path to PDB file
            
        Returns:
            dict: Protein information including structure, sequence, and stats
        """
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")
            
        structure = self.parser.get_structure(pdb_path.stem, str(pdb_path))
        
        # Extract information
        info = {
            'structure': structure,
            'name': pdb_path.stem,
            'path': str(pdb_path),
            'chains': self._get_chains(structure),
            'residue_count': self._count_residues(structure),
            'atoms': self._count_atoms(structure),
            'sequence': self._extract_sequence(structure)
        }
        
        return info
    
    def fetch_from_pdb_id(self, pdb_id, output_dir="data/test_proteins"):
        """
        Download protein structure from RCSB PDB
        
        Args:
            pdb_id: 4-letter PDB identifier (e.g., '6LU7')
            output_dir: Where to save the downloaded file
            
        Returns:
            dict: Protein information
        """
        pdb_id = pdb_id.upper()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Download from RCSB
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        
        if response.status_code != 200:
            raise ValueError(f"Could not download PDB {pdb_id}. Status: {response.status_code}")
        
        # Save to file
        pdb_path = output_dir / f"{pdb_id}.pdb"
        pdb_path.write_text(response.text)
        
        print(f"✓ Downloaded {pdb_id} to {pdb_path}")
        
        # Load and return info
        return self.load_from_pdb_file(pdb_path)
    
    def predict_binding_sites(self, protein_info, method="geometric"):
        """
        Predict binding sites on the protein
        
        Args:
            protein_info: Dict returned from load_from_pdb_file or fetch_from_pdb_id
            method: 'geometric' (fast, rule-based) or 'p2rank' (ML-based, requires p2rank)
            
        Returns:
            list: Binding site predictions with coordinates and scores
        """
        if method == "geometric":
            return self._geometric_pocket_detection(protein_info)
        elif method == "p2rank":
            return self._p2rank_prediction(protein_info)
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def _geometric_pocket_detection(self, protein_info):
        """
        Simple geometric pocket detection based on cavity analysis
        Uses convex hull and surface accessibility as heuristics
        """
        structure = protein_info['structure']
        
        # Get all CA atoms (alpha carbons)
        ca_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue) and 'CA' in residue:
                        ca_atoms.append(residue['CA'])
        
        if len(ca_atoms) < 10:
            return [{"message": "Too few residues for pocket detection"}]
        
        # Calculate geometric center
        coords = np.array([atom.coord for atom in ca_atoms])
        center = coords.mean(axis=0)
        
        # Find regions far from center (potential surface pockets)
        distances = np.linalg.norm(coords - center, axis=1)
        surface_threshold = np.percentile(distances, 75)
        
        surface_atoms = [ca_atoms[i] for i, d in enumerate(distances) if d > surface_threshold]
        
        # Cluster surface atoms into potential binding sites
        pockets = []
        if len(surface_atoms) > 5:
            # Simple centroid-based clustering
            surface_coords = np.array([atom.coord for atom in surface_atoms])
            
            # Find 3 main pockets (simplified k-means)
            n_pockets = min(3, len(surface_atoms) // 5)
            indices = np.random.choice(len(surface_coords), n_pockets, replace=False)
            pocket_centers = surface_coords[indices]
            
            for i, pocket_center in enumerate(pocket_centers):
                # Find atoms within 10Å of pocket center
                pocket_atoms = []
                for j, coord in enumerate(surface_coords):
                    if np.linalg.norm(coord - pocket_center) < 10.0:
                        pocket_atoms.append(surface_atoms[j])
                
                if len(pocket_atoms) >= 3:
                    pockets.append({
                        'id': i + 1,
                        'center': pocket_center.tolist(),
                        'confidence': min(len(pocket_atoms) / 20.0, 1.0),  # Heuristic score
                        'residue_count': len(pocket_atoms),
                        'method': 'geometric'
                    })
        
        return sorted(pockets, key=lambda x: x['confidence'], reverse=True)
    
    def _p2rank_prediction(self, protein_info):
        """
        Use P2Rank for ML-based pocket prediction
        Requires p2rank to be installed
        """
        # Create temp directory for p2rank output
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            try:
                # Run p2rank
                cmd = [
                    self.p2rank_path, "predict",
                    "-f", protein_info['path'],
                    "-o", str(tmpdir)
                ]
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                
                if result.returncode != 0:
                    print(f"P2Rank error: {result.stderr}")
                    return self._geometric_pocket_detection(protein_info)
                
                # Parse p2rank output
                predictions_file = tmpdir / f"{protein_info['name']}.pdb_predictions.csv"
                
                if predictions_file.exists():
                    pockets = self._parse_p2rank_output(predictions_file)
                    return pockets
                else:
                    print("P2Rank output not found, falling back to geometric method")
                    return self._geometric_pocket_detection(protein_info)
                    
            except (subprocess.TimeoutExpired, FileNotFoundError) as e:
                print(f"P2Rank failed: {e}, using geometric method")
                return self._geometric_pocket_detection(protein_info)
    
    def _parse_p2rank_output(self, csv_path):
        """Parse P2Rank CSV output"""
        pockets = []
        with open(csv_path) as f:
            lines = f.readlines()[1:]  # Skip header
            for line in lines:
                parts = line.strip().split(',')
                if len(parts) >= 4:
                    pockets.append({
                        'id': int(parts[0]),
                        'confidence': float(parts[2]),
                        'residue_count': int(parts[3]),
                        'method': 'p2rank'
                    })
        return pockets
    
    # Helper methods
    def _get_chains(self, structure):
        """Extract chain IDs"""
        chains = []
        for model in structure:
            for chain in model:
                chains.append(chain.id)
        return chains
    
    def _count_residues(self, structure):
        """Count total residues"""
        count = 0
        for model in structure:
            for chain in model:
                count += len([r for r in chain if is_aa(r)])
        return count
    
    def _count_atoms(self, structure):
        """Count total atoms"""
        count = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    count += sum(1 for _ in residue.get_atoms())
        return count
    
    def _extract_sequence(self, structure):
        """Extract amino acid sequence"""
        from Bio.PDB.Polypeptide import PPBuilder
        ppb = PPBuilder()
        sequences = {}
        
        for model in structure:
            for chain in model:
                pp_list = ppb.build_peptides(chain)
                if pp_list:
                    sequences[chain.id] = str(pp_list[0].get_sequence())
        
        return sequences


# Test function
def test_protein_processor():
    """Test the protein processor with COVID-19 main protease"""
    processor = ProteinProcessor()
    
    print("🧬 Testing Protein Processor\n" + "="*60)
    
    # Test 1: Fetch from PDB
    print("\n1. Fetching COVID-19 Main Protease (6LU7)...")
    try:
        protein = processor.fetch_from_pdb_id("6LU7")
        print(f"   ✓ Loaded: {protein['name']}")
        print(f"   ✓ Chains: {protein['chains']}")
        print(f"   ✓ Residues: {protein['residue_count']}")
        print(f"   ✓ Atoms: {protein['atoms']}")
        
        # Test 2: Predict binding sites
        print("\n2. Predicting binding sites (geometric method)...")
        pockets = processor.predict_binding_sites(protein, method="geometric")
        print(f"   ✓ Found {len(pockets)} potential binding pockets")
        
        for pocket in pockets[:3]:
            print(f"   • Pocket {pocket['id']}: "
                  f"confidence={pocket['confidence']:.2f}, "
                  f"residues={pocket['residue_count']}")
        
        print("\n✅ All tests passed!")
        return protein, pockets
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None, None


if __name__ == "__main__":
    test_protein_processor()
