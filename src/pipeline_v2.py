"""
PharmaLens Pipeline - End-to-End Drug Discovery
Combines protein processing, molecule generation, docking, and literature search
"""

from protein_processor import ProteinProcessor
from docking_engine import DockingEngine
from molecule_generator import MoleculeGenerator
from literature_search import LiteratureSearch
import time


def run_drug_discovery_pipeline(pdb_id="6LU7", n_molecules=10, top_n=5, 
                                include_literature=True):
    """
    Complete drug discovery pipeline with literature search
    
    Args:
        pdb_id: PDB identifier for target protein
        n_molecules: Number of molecules to generate
        top_n: Number of top candidates to return
        include_literature: Whether to search PubMed for relevant papers
        
    Returns:
        list: Top docking results with literature citations
    """
    print("🧬 PharmaLens Drug Discovery Pipeline")
    print("="*70)
    
    # Step 1: Load and analyze protein
    print(f"\n[1/5] Loading target protein ({pdb_id})...")
    processor = ProteinProcessor()
    protein = processor.fetch_from_pdb_id(pdb_id)
    print(f"      ✓ Loaded {protein['name']}")
    print(f"      ✓ Chains: {protein['chains']}")
    print(f"      ✓ Residues: {protein['residue_count']}")
    
    # Step 2: Find binding sites
    print(f"\n[2/5] Predicting binding sites...")
    pockets = processor.predict_binding_sites(protein, method="geometric")
    best_pocket = pockets[0]
    print(f"      ✓ Found {len(pockets)} pockets")
    print(f"      ✓ Best pocket: {best_pocket['id']} "
          f"(confidence: {best_pocket['confidence']:.2f})")
    print(f"      ✓ Center: [{best_pocket['center'][0]:.1f}, "
          f"{best_pocket['center'][1]:.1f}, {best_pocket['center'][2]:.1f}]")
    
    # Step 3: Generate molecules
    print(f"\n[3/5] Generating {n_molecules} drug candidates...")
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
    print(f"\n[4/5] Docking molecules to binding pocket...")
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
                    'druglike': props['molecular_weight'] <= 500 and props['logP'] <= 5,
                    'papers': []  # Will be populated in next step
                })
        
        except Exception as e:
            # Silently continue - some molecules may fail
            continue
    
    elapsed_time = time.time() - start_time
    
    # Sort by binding score (lower is better)
    results.sort(key=lambda x: x['binding_score'])
    
    # Step 5: Literature search for top candidates
    if include_literature and results:
        print(f"\n[5/5] Searching literature for top candidates...")
        searcher = LiteratureSearch()
        
        # Get target name
        target_names = {
            "6LU7": "COVID-19 main protease",
            "1HSG": "HIV-1 protease",
            "4AGN": "EGFR kinase"
        }
        target_name = target_names.get(pdb_id, "protein target")
        
        for i, result in enumerate(results[:top_n], 1):
            try:
                print(f"      • Finding papers for candidate {i}...")
                
                # Search for papers about similar molecules and the target
                papers = searcher.search_by_target(
                    target_name=target_name,
                    max_results=3
                )
                
                result['papers'] = papers
                
            except Exception as e:
                print(f"      ⚠️  Literature search failed for candidate {i}")
                continue
        
        print(f"      ✓ Found literature for top {min(top_n, len(results))} candidates")
    else:
        print(f"\n[5/5] Skipping literature search")
    
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
        
        # Show literature
        if result['papers']:
            print(f"   📚 Relevant Papers:")
            for paper in result['papers'][:2]:  # Show top 2
                print(f"      • {paper['title'][:60]}...")
                print(f"        {paper['authors']} ({paper['year']})")
                print(f"        {paper['url']}")
        
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
    
    # Run the full pipeline with literature search
    results = run_drug_discovery_pipeline(
        pdb_id="6LU7",
        n_molecules=10,
        top_n=5,
        include_literature=True
    )
    
    # Benchmark against known drug
    compare_to_known_drug(pdb_id="6LU7", drug_name="Nirmatrelvir")
    
    print("\n" + "="*70)
    print("🎉 Drug Discovery Complete!")
    print("="*70)
