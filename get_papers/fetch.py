from typing import List
from Bio import Entrez
import time

# Set your real email here (required by NCBI)
Entrez.email = "your_email@example.com"

def search_pubmed(query: str, max_results: int = 50) -> List[str]:
    """Search PubMed for a query and return a list of paper IDs."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    results = Entrez.read(handle)
    handle.close()
    return results["IdList"]

def fetch_details(paper_ids: List[str]) -> str:
    """Fetch details for a list of PubMed IDs in MEDLINE text format."""
    time.sleep(1)  # Avoid API rate limiting
    handle = Entrez.efetch(
        db="pubmed",
        id=",".join(paper_ids),
        rettype="medline",
        retmode="text"
    )
    data = handle.read()
    handle.close()
    return data
