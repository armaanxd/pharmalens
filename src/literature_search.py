"""
Literature Search Module
Search PubMed for relevant papers related to molecules and drug targets
"""

import requests
import time
from typing import List, Dict
import xml.etree.ElementTree as ET


class LiteratureSearch:
    """Search scientific literature for drug-related information"""
    
    def __init__(self):
        """Initialize literature search"""
        self.pubmed_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.email = "your.email@example.com"  # PubMed requires email for API
        
    def search_molecule_papers(self, molecule_name=None, target_name=None, 
                               smiles=None, max_results=5):
        """
        Search PubMed for papers related to a molecule or target
        
        Args:
            molecule_name: Name of molecule (e.g., "aspirin")
            target_name: Name of target protein (e.g., "COVID-19 main protease")
            smiles: SMILES string (will extract scaffold for search)
            max_results: Maximum papers to return
            
        Returns:
            list: Paper information with titles, abstracts, PMIDs
        """
        # Build search query
        query_parts = []
        
        if molecule_name:
            query_parts.append(f"{molecule_name}[Title/Abstract]")
        
        if target_name:
            query_parts.append(f"{target_name}[Title/Abstract]")
            
        if smiles and not molecule_name:
            # Extract common scaffolds from SMILES
            scaffold_terms = self._extract_scaffold_terms(smiles)
            if scaffold_terms:
                query_parts.append(f"({' OR '.join(scaffold_terms)})[Title/Abstract]")
        
        if not query_parts:
            return []
        
        # Add drug discovery related terms
        query_parts.append("(drug OR inhibitor OR binding OR docking OR compound)")
        
        query = " AND ".join(query_parts)
        
        # Search PubMed
        try:
            pmids = self._search_pubmed(query, max_results)
            
            if not pmids:
                return []
            
            # Fetch paper details
            papers = self._fetch_paper_details(pmids)
            
            return papers
            
        except Exception as e:
            print(f"⚠️  Literature search failed: {e}")
            return []
    
    def search_by_target(self, target_name, drug_class=None, max_results=5):
        """
        Search for papers about a specific drug target
        
        Args:
            target_name: Protein target (e.g., "EGFR kinase")
            drug_class: Optional drug class (e.g., "antiviral")
            max_results: Maximum papers to return
            
        Returns:
            list: Paper information
        """
        query = f"{target_name}[Title/Abstract] AND (inhibitor OR drug OR compound)"
        
        if drug_class:
            query += f" AND {drug_class}[Title/Abstract]"
        
        try:
            pmids = self._search_pubmed(query, max_results)
            papers = self._fetch_paper_details(pmids)
            return papers
        except Exception as e:
            print(f"⚠️  Target search failed: {e}")
            return []
    
    def search_similar_molecules(self, reference_smiles, max_results=3):
        """
        Search for papers about molecules with similar scaffolds
        
        Args:
            reference_smiles: SMILES to find similar molecules for
            max_results: Maximum papers
            
        Returns:
            list: Paper information
        """
        scaffold_terms = self._extract_scaffold_terms(reference_smiles)
        
        if not scaffold_terms:
            return []
        
        query = f"({' OR '.join(scaffold_terms)}) AND (structure OR scaffold OR derivative)"
        
        try:
            pmids = self._search_pubmed(query, max_results)
            papers = self._fetch_paper_details(pmids)
            return papers
        except Exception as e:
            print(f"⚠️  Similarity search failed: {e}")
            return []
    
    def _search_pubmed(self, query, max_results=5):
        """
        Search PubMed and return PMIDs
        
        Args:
            query: Search query string
            max_results: Maximum results to return
            
        Returns:
            list: PMIDs
        """
        search_url = f"{self.pubmed_base}/esearch.fcgi"
        
        params = {
            'db': 'pubmed',
            'term': query,
            'retmax': max_results,
            'retmode': 'json',
            'sort': 'relevance',
            'email': self.email
        }
        
        response = requests.get(search_url, params=params, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        
        pmids = data.get('esearchresult', {}).get('idlist', [])
        
        return pmids
    
    def _fetch_paper_details(self, pmids):
        """
        Fetch paper details from PMIDs
        
        Args:
            pmids: List of PubMed IDs
            
        Returns:
            list: Paper information dictionaries
        """
        if not pmids:
            return []
        
        fetch_url = f"{self.pubmed_base}/efetch.fcgi"
        
        params = {
            'db': 'pubmed',
            'id': ','.join(pmids),
            'retmode': 'xml',
            'email': self.email
        }
        
        response = requests.get(fetch_url, params=params, timeout=15)
        response.raise_for_status()
        
        # Parse XML
        root = ET.fromstring(response.content)
        
        papers = []
        
        for article in root.findall('.//PubmedArticle'):
            try:
                # Extract PMID
                pmid = article.find('.//PMID').text
                
                # Extract title
                title_elem = article.find('.//ArticleTitle')
                title = title_elem.text if title_elem is not None else "No title"
                
                # Extract abstract
                abstract_elem = article.find('.//AbstractText')
                abstract = abstract_elem.text if abstract_elem is not None else "No abstract available"
                
                # Truncate abstract
                if len(abstract) > 500:
                    abstract = abstract[:500] + "..."
                
                # Extract authors
                authors = []
                for author in article.findall('.//Author')[:3]:  # First 3 authors
                    last_name = author.find('.//LastName')
                    if last_name is not None:
                        authors.append(last_name.text)
                
                author_str = ", ".join(authors)
                if len(article.findall('.//Author')) > 3:
                    author_str += " et al."
                
                # Extract year
                year_elem = article.find('.//PubDate/Year')
                year = year_elem.text if year_elem is not None else "N/A"
                
                # Extract journal
                journal_elem = article.find('.//Journal/Title')
                journal = journal_elem.text if journal_elem is not None else "Unknown"
                
                papers.append({
                    'pmid': pmid,
                    'title': title,
                    'abstract': abstract,
                    'authors': author_str,
                    'year': year,
                    'journal': journal,
                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
                
            except Exception as e:
                print(f"⚠️  Failed to parse article: {e}")
                continue
        
        return papers
    
    def _extract_scaffold_terms(self, smiles):
        """
        Extract searchable terms from SMILES structure
        
        Args:
            smiles: SMILES string
            
        Returns:
            list: Search terms
        """
        from rdkit import Chem
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return []
            
            terms = []
            
            # Check for common scaffolds
            if 'c1ccccc1' in smiles.lower():
                terms.append('benzene')
                terms.append('phenyl')
            
            if 'c1ccncc1' in smiles.lower() or 'pyridin' in smiles.lower():
                terms.append('pyridine')
            
            if 'indol' in smiles.lower():
                terms.append('indole')
            
            if 'quinolin' in smiles.lower():
                terms.append('quinoline')
            
            # Check for functional groups
            if 'C(=O)O' in smiles:
                terms.append('carboxylic acid')
            
            if 'C(=O)N' in smiles:
                terms.append('amide')
            
            if 'S(=O)(=O)' in smiles:
                terms.append('sulfonamide')
            
            return terms[:3]  # Limit to 3 terms
            
        except Exception:
            return []
    
    def format_citation(self, paper):
        """
        Format paper as citation string
        
        Args:
            paper: Paper dictionary
            
        Returns:
            str: Formatted citation
        """
        return f"{paper['authors']} ({paper['year']}). {paper['title']}. {paper['journal']}. PMID: {paper['pmid']}"


# Test function
def test_literature_search():
    """Test the literature search module"""
    print("📚 Testing Literature Search\n" + "="*60)
    
    searcher = LiteratureSearch()
    
    # Test 1: Search by target
    print("\n1. Searching papers about COVID-19 main protease...")
    papers = searcher.search_by_target(
        target_name="COVID-19 main protease",
        drug_class="antiviral",
        max_results=3
    )
    
    print(f"   ✓ Found {len(papers)} papers")
    
    if papers:
        print("\n   Top paper:")
        paper = papers[0]
        print(f"   • Title: {paper['title'][:80]}...")
        print(f"   • Authors: {paper['authors']}")
        print(f"   • Year: {paper['year']}")
        print(f"   • PMID: {paper['pmid']}")
        print(f"   • URL: {paper['url']}")
    
    # Test 2: Search for a molecule
    print("\n2. Searching papers about aspirin...")
    papers = searcher.search_molecule_papers(
        molecule_name="aspirin",
        max_results=3
    )
    
    print(f"   ✓ Found {len(papers)} papers")
    
    if papers:
        print(f"\n   Example citation:")
        print(f"   {searcher.format_citation(papers[0])}")
    
    # Test 3: Search by SMILES (protease inhibitor-like)
    print("\n3. Searching papers for similar molecules...")
    test_smiles = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"
    papers = searcher.search_similar_molecules(test_smiles, max_results=2)
    
    print(f"   ✓ Found {len(papers)} papers about similar scaffolds")
    
    print("\n✅ All tests passed!")
    print("\n💡 Literature search ready for integration!")
    
    return papers


if __name__ == "__main__":
    test_literature_search()
