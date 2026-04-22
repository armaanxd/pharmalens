"""
PharmaLens Pipeline - End-to-End Drug Discovery
Combines protein processing, molecule generation, and docking
"""

from protein_processor import ProteinProcessor
from docking_engine import DockingEngine
from molecule_generator import MoleculeGenerator
import time


def run_drug_discovery_pipeline(pdb_id="6LU7", n_molecules=10, top_n=5):
    """
    Complete drug discovery pipeline
    
    Args:
        pdb_id: PDB identifier for target protein
        n_molecules: Number of molecules to generate
        top_n: Number of top candidates to return
        
    Returns:
        list: Top docking results
    """
    print("🧬 PharmaLens Drug Discovery Pipeline")
    print("="*70)
    
    # Step 1: Load and analyze protein
    print(f"\n[1/4] Loading target protein ({pdb_id})...")
    processor = ProteinProcessor()
    protein = processor.fetch_from_pdb_id(pdb_id)
    print(f"      ✓ Loaded {protein['name']}")
    print(f"      ✓ Chains: {protein['chains']}")
    print(f"      ✓ Residues: {protein['residue_count']}")
    
    # Step 2: Find binding sites
    print(f"\n[2/4] Predicting binding sites...")
    pockets = processor.predict_binding_sites(protein, method="geometric")
    best_pocket = pockets[0]
    print(f"      ✓ Found {len(pockets)} pockets")
    print(f"      ✓ Best pocket: {best_pocket['id']} "
          f"(confidence: {best_pocket['confidence']:.2f})")
    print(f"      ✓ Center: [{best_pocket['center'][0]:.1f}, "
          f"{best_pocket['center'][1]:.1f}, {best_pocket['center'][2]:.1f}]")
    
    # Step 3: Generate molecules
    print(f"\n[3/4] Generating {n_molecules} drug candidates...")
    generator = MoleculeGenerator()
    
    # Generate both random and template-based molecules
    antiviral_mols = generator.generate_based_on_known_drugs(
        drug_class="antiviral",
        n_molecules=n_molecules // 2
    )
    random_mols = generator.generate_random_molecules(
        n_molecules=n_molecules // 2
    )
    
    all_molecules = antiviral_mols + random_mols
    
    # Filter for druglikeness
    filtered_molecules = generator.filter_by_properties(
        all_molecules,
        min_mw=200,
        max_mw=500,
        min_logp=-1,
        max_logp=5
    )
    
    print(f"      ✓ Generated {len(all_molecules)} molecules")
    print(f"      ✓ {len(filtered_molecules)} pass druglikeness filters")
    
    # Step 4: Dock molecules
    print(f"\n[4/4] Docking molecules to binding pocket...")
    print(f"      (This may take a few minutes...)")
    
    docker = DockingEngine()
    results = []
    
    start_time = time.time()
    
    for i, smiles in enumerate(filtered_molecules, 1):
        try:
            # Show progress
            if i % 3 == 0 or i == 1:
                print(f"      • Docking molecule {i}/{len(filtered_molecules)}...")
            
            # Dock molecule
            docking_result = docker.dock_molecule(
                protein,
                smiles,
                best_pocket['center'],
                exhaustiveness=4  # Lower for speed
            )
            
            if docking_result.get('success'):
                # Get properties
                props = generator.get_molecule_properties(smiles)
                
                # Combine results
                results.append({
                    'smiles': smiles,
                    'binding_score': docking_result['best_score'],
                    'molecular_weight': props['molecular_weight'],
                    'logP': props['logP'],
                    'h_donors': props['h_bond_donors'],
                    'h_acceptors': props['h_bond_acceptors'],
                    'druglike': props['molecular_weight'] <= 500 and props['logP'] <= 5
                })
        
        except Exception as e:
            print(f"      ⚠️  Molecule {i} failed: {e}")
            continue
    
    elapsed_time = time.time() - start_time
    
    # Sort by binding score (lower is better)
    results.sort(key=lambda x: x['binding_score'])
    
    # Display results
    print(f"\n{'='*70}")
    print(f"✅ Pipeline Complete! ({elapsed_time:.1f} seconds)")
    print(f"{'='*70}")
    
    print(f"\n🏆 Top {min(top_n, len(results))} Drug Candidates:\n")
    
    for i, result in enumerate(results[:top_n], 1):
        print(f"{i}. Binding Score: {result['binding_score']:.2f} kcal/mol")
        print(f"   Molecular Weight: {result['molecular_weight']}")
        print(f"   LogP: {result['logP']}")
        print(f"   H-Bond Donors/Acceptors: {result['h_donors']}/{result['h_acceptors']}")
        print(f"   SMILES: {result['smiles'][:60]}...")
        print(f"   Druglike: {'✓ Yes' if result['druglike'] else '✗ No'}")
        print()
    
    # Save results
    save_results(results, pdb_id)
    
    return results


def save_results(results, pdb_id):
    """Save pipeline results to file"""
    from pathlib import Path
    import json
    
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)
    
    # Save as JSON
    output_file = output_dir / f"{pdb_id}_drug_candidates.json"
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"💾 Results saved to: {output_file}")
    
    # Also save top SMILES to text file
    smiles_file = output_dir / f"{pdb_id}_top_molecules.txt"
    with open(smiles_file, 'w') as f:
        for i, result in enumerate(results[:10], 1):
            f.write(f"{i}. Score: {result['binding_score']:.2f} | "
                   f"{result['smiles']}\n")
    
    print(f"💾 Top SMILES saved to: {smiles_file}")


def compare_to_known_drug(pdb_id="6LU7", drug_name="Nirmatrelvir"):
    """
    Compare pipeline results to a known drug
    """
    print(f"\n{'='*70}")
    print(f"🔬 Benchmark: Comparing to Known Drug ({drug_name})")
    print(f"{'='*70}\n")
    
    # Known drugs and their SMILES
    known_drugs = {
        "Nirmatrelvir": "CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C",
        "Remdesivir": "CCC(CC)COC(=O)C(C)NP(=O)(OCC1C(C(C(O1)N2C=CC(=NC2=O)N)O)O)OC3=CC=CC=C3",
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O"
    }
    
    if drug_name not in known_drugs:
        print(f"⚠️  Drug {drug_name} not in database")
        return
    
    drug_smiles = known_drugs[drug_name]
    
    # Load protein
    processor = ProteinProcessor()
    protein = processor.fetch_from_pdb_id(pdb_id)
    pockets = processor.predict_binding_sites(protein)
    
    # Dock known drug
    docker = DockingEngine()
    generator = MoleculeGenerator()
    
    print(f"Docking {drug_name}...")
    result = docker.dock_molecule(
        protein,
        drug_smiles,
        pockets[0]['center']
    )
    
    if result.get('success'):
        props = generator.get_molecule_properties(drug_smiles)
        
        print(f"\n{drug_name} Results:")
        print(f"  Binding Score: {result['best_score']:.2f} kcal/mol")
        print(f"  Molecular Weight: {props['molecular_weight']}")
        print(f"  LogP: {props['logP']}")
        print(f"  SMILES: {drug_smiles[:60]}...")
        
        print(f"\n💡 Use this as a benchmark for your AI-generated molecules!")
        print(f"   Goal: Find molecules with binding scores < {result['best_score']:.2f}")


if __name__ == "__main__":
    print("Starting PharmaLens Drug Discovery Pipeline...\n")
    
    # Run the full pipeline
    results = run_drug_discovery_pipeline(
        pdb_id="6LU7",
        n_molecules=10,
        top_n=5
    )
    
    # Benchmark against known drug
    compare_to_known_drug(pdb_id="6LU7", drug_name="Nirmatrelvir")
    
    print("\n" + "="*70)
    print("🎉 Drug Discovery Complete!")
    print("="*70)
