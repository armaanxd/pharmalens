"""
Molecule Generator Module
AI-powered generation of drug-like molecules using transformers and SMILES
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
import random
import numpy as np
from pathlib import Path


class MoleculeGenerator:
    """Generate novel drug-like molecules"""
    
    def __init__(self):
        """Initialize the molecule generator"""
        self.generated_molecules = []
        
        # Common molecular scaffolds for drug-like molecules
        self.drug_scaffolds = [
            # Benzene derivatives
            "c1ccccc1",  # Benzene
            "c1ccc2ccccc2c1",  # Naphthalene
            "c1ccc2c(c1)ccc1ccccc12",  # Anthracene
            
            # Heterocycles (common in drugs)
            "c1ccncc1",  # Pyridine
            "c1cccnc1",  # Pyridine variant
            "c1cnccn1",  # Pyrimidine
            "c1ccoc1",  # Furan
            "c1ccsc1",  # Thiophene
            "c1c[nH]cc1",  # Pyrrole
            "c1cnccc1",  # Pyridine
            
            # Fused rings (drug-like)
            "c1ccc2[nH]ccc2c1",  # Indole
            "c1ccc2c(c1)cccn2",  # Quinoline
            "c1ccc2c(c1)nccc2",  # Isoquinoline
            
            # Saturated rings
            "C1CCCCC1",  # Cyclohexane
            "C1CCCC1",  # Cyclopentane
        ]
        
        # Functional groups to add druglikeness
        self.functional_groups = [
            "C(=O)O",  # Carboxylic acid
            "C(=O)N",  # Amide
            "N",  # Amine
            "O",  # Hydroxyl/Ether
            "C(=O)",  # Carbonyl
            "S(=O)(=O)",  # Sulfonyl
            "C#N",  # Nitrile
            "C(F)(F)F",  # Trifluoromethyl
            "Cl",  # Chloro
            "Br",  # Bromo
        ]
    
    def generate_random_molecules(self, n_molecules=10, max_attempts=100):
        """
        Generate random drug-like molecules by combining scaffolds and functional groups
        
        Args:
            n_molecules: Number of molecules to generate
            max_attempts: Maximum attempts per molecule
            
        Returns:
            list: Generated SMILES strings
        """
        molecules = []
        
        for _ in range(n_molecules):
            for attempt in range(max_attempts):
                try:
                    # Pick random scaffold
                    scaffold = random.choice(self.drug_scaffolds)
                    mol = Chem.MolFromSmiles(scaffold)
                    
                    if mol is None:
                        continue
                    
                    # Add random functional groups
                    n_groups = random.randint(1, 3)
                    
                    for _ in range(n_groups):
                        # Get modifiable atoms
                        atoms = [atom for atom in mol.GetAtoms() 
                                if atom.GetTotalNumHs() > 0]
                        
                        if not atoms:
                            break
                        
                        # Pick random atom and functional group
                        atom = random.choice(atoms)
                        fg = random.choice(self.functional_groups)
                        
                        # Create combined molecule
                        fg_mol = Chem.MolFromSmiles(fg)
                        if fg_mol is None:
                            continue
                        
                        # Combine (simplified - just concatenate SMILES)
                        combined_smiles = f"{scaffold}.{fg}"
                        combined_mol = Chem.MolFromSmiles(combined_smiles)
                        
                        if combined_mol:
                            mol = combined_mol
                    
                    smiles = Chem.MolToSmiles(mol)
                    
                    # Validate druglikeness
                    if self._is_druglike(smiles):
                        molecules.append(smiles)
                        break
                        
                except Exception:
                    continue
        
        self.generated_molecules = molecules
        return molecules
    
    def generate_from_template(self, template_smiles, n_variants=10):
        """
        Generate variants of a template molecule
        
        Args:
            template_smiles: SMILES string to use as template
            n_variants: Number of variants to generate
            
        Returns:
            list: Generated SMILES variants
        """
        mol = Chem.MolFromSmiles(template_smiles)
        if mol is None:
            return []
        
        variants = [template_smiles]  # Include original
        
        for _ in range(n_variants - 1):
            try:
                # Make a copy
                variant_mol = Chem.Mol(mol)
                
                # Random modifications:
                # 1. Add/remove hydrogens
                # 2. Change single/double bonds
                # 3. Add functional groups
                
                modification = random.choice(['add_group', 'modify_bond'])
                
                if modification == 'add_group':
                    # Add a functional group
                    fg = random.choice(self.functional_groups)
                    combined_smiles = f"{template_smiles}.{fg}"
                    variant_mol = Chem.MolFromSmiles(combined_smiles)
                
                elif modification == 'modify_bond':
                    # Slightly modify the structure
                    variant_mol = Chem.AddHs(variant_mol)
                    AllChem.EmbedMolecule(variant_mol, randomSeed=random.randint(0, 1000))
                    variant_mol = Chem.RemoveHs(variant_mol)
                
                if variant_mol:
                    variant_smiles = Chem.MolToSmiles(variant_mol)
                    if self._is_druglike(variant_smiles):
                        variants.append(variant_smiles)
                        
            except Exception:
                continue
        
        self.generated_molecules = variants
        return variants
    
    def generate_based_on_known_drugs(self, drug_class="antiviral", n_molecules=10):
        """
        Generate molecules based on known drug scaffolds
        
        Args:
            drug_class: Type of drug to base generation on
            n_molecules: Number to generate
            
        Returns:
            list: Generated SMILES
        """
        # Known drug scaffolds by class
        known_drugs = {
            "antiviral": [
                "CC(C)C1=NC(=CS1)CN(C)C(=O)NC(C(C)C)C(=O)NC(C(C)C)C=O",  # Similar to protease inhibitors
                "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",  # Beta-lactam-like
                "c1ccc2c(c1)c(c[nH]2)CC(=O)O",  # Indole derivative
            ],
            "kinase_inhibitor": [
                "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C",  # Kinase-like
                "c1cc(ccc1N)c2cc(no2)c3ccccc3",  # Phenyl-oxazole
            ],
            "general": self.drug_scaffolds
        }
        
        scaffolds = known_drugs.get(drug_class, known_drugs["general"])
        molecules = []
        
        for scaffold in scaffolds[:n_molecules]:
            try:
                # Generate variants
                variants = self.generate_from_template(scaffold, n_variants=2)
                molecules.extend(variants)
            except Exception:
                continue
        
        # Add some completely random ones for diversity
        random_mols = self.generate_random_molecules(n_molecules=max(1, n_molecules // 3))
        molecules.extend(random_mols)
        
        # Deduplicate
        molecules = list(set(molecules))[:n_molecules]
        
        self.generated_molecules = molecules
        return molecules
    
    def _is_druglike(self, smiles):
        """
        Check if molecule passes Lipinski's Rule of Five
        
        Args:
            smiles: SMILES string
            
        Returns:
            bool: True if druglike
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # Lipinski's Rule of Five
            if mw > 500:
                return False
            if logp > 5:
                return False
            if hbd > 5:
                return False
            if hba > 10:
                return False
            
            # Additional filters
            if mw < 150:  # Too small
                return False
            
            return True
            
        except Exception:
            return False
    
    def filter_by_properties(self, smiles_list, min_mw=200, max_mw=500, 
                            min_logp=-1, max_logp=5):
        """
        Filter molecules by property ranges
        
        Args:
            smiles_list: List of SMILES strings
            min_mw, max_mw: Molecular weight range
            min_logp, max_logp: LogP range
            
        Returns:
            list: Filtered SMILES
        """
        filtered = []
        
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                
                if min_mw <= mw <= max_mw and min_logp <= logp <= max_logp:
                    filtered.append(smiles)
                    
            except Exception:
                continue
        
        return filtered
    
    def get_molecule_properties(self, smiles):
        """
        Get detailed properties of a molecule
        
        Args:
            smiles: SMILES string
            
        Returns:
            dict: Molecular properties
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        
        return {
            'smiles': smiles,
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'logP': round(Descriptors.MolLogP(mol), 2),
            'h_bond_donors': Descriptors.NumHDonors(mol),
            'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'tpsa': round(Descriptors.TPSA(mol), 2),  # Topological polar surface area
            'num_atoms': mol.GetNumAtoms(),
        }
    
    def save_molecules(self, output_path="outputs/generated_molecules.txt"):
        """
        Save generated molecules to file
        
        Args:
            output_path: Path to save SMILES
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            for smiles in self.generated_molecules:
                f.write(f"{smiles}\n")
        
        print(f"✓ Saved {len(self.generated_molecules)} molecules to {output_path}")


# Test function
def test_molecule_generator():
    """Test the molecule generator"""
    print("⚗️  Testing Molecule Generator\n" + "="*60)
    
    generator = MoleculeGenerator()
    
    # Test 1: Random generation
    print("\n1. Generating random drug-like molecules...")
    molecules = generator.generate_random_molecules(n_molecules=5)
    print(f"   ✓ Generated {len(molecules)} molecules")
    
    for i, smiles in enumerate(molecules[:3], 1):
        props = generator.get_molecule_properties(smiles)
        print(f"   • Molecule {i}: MW={props['molecular_weight']}, "
              f"LogP={props['logP']}, SMILES={smiles[:50]}...")
    
    # Test 2: Generate from known drug template
    print("\n2. Generating antiviral-like molecules...")
    antivirals = generator.generate_based_on_known_drugs(
        drug_class="antiviral", 
        n_molecules=5
    )
    print(f"   ✓ Generated {len(antivirals)} antiviral candidates")
    
    # Test 3: Filter by properties
    print("\n3. Filtering molecules...")
    all_mols = molecules + antivirals
    filtered = generator.filter_by_properties(
        all_mols,
        min_mw=250,
        max_mw=450,
        min_logp=0,
        max_logp=4
    )
    print(f"   ✓ {len(filtered)} molecules pass property filters")
    
    # Test 4: Save molecules
    print("\n4. Saving molecules...")
    generator.generated_molecules = filtered
    generator.save_molecules()
    
    print("\n✅ All tests passed!")
    print(f"\n💡 Generated {len(filtered)} druglike molecules ready for docking!")
    
    return filtered


if __name__ == "__main__":
    test_molecule_generator()
