import argparse
import csv
import sys
from Bio import Entrez
import re

# Set your email for NCBI Entrez API
Entrez.email = "your_email@example.com"  # Change this to your actual email

# List of keywords to detect company affiliations
COMPANY_KEYWORDS = ["pharma", "biotech", "therapeutics", "laboratories", "inc", "ltd", "gmbh", "corporation"]

def fetch_pubmed_papers(query, debug=False):
    """Fetches PubMed papers based on the given query."""
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
        record = Entrez.read(handle)
        handle.close()
        pubmed_ids = record.get("IdList", [])
        
        papers = []
        for pubmed_id in pubmed_ids:
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="medline", retmode="text")
            article = handle.read()
            handle.close()
            
            title_match = re.search(r'TI  - (.+)', article)
            title = title_match.group(1) if title_match else "Unknown Title"
            
            date_match = re.search(r'DP  - (\d{4})', article)
            publication_date = date_match.group(1) if date_match else "Unknown Date"
            
            authors = re.findall(r'AU  - (.+)', article)
            affiliations = re.findall(r'AD  - (.+)', article)
            
            non_academic_authors = []
            company_affiliations = []
            
            for author, affil in zip(authors, affiliations):
                if any(keyword in affil.lower() for keyword in COMPANY_KEYWORDS):
                    non_academic_authors.append(author)
                    company_affiliations.append(affil)
            
            email_match = re.search(r'FAU - .*?\n(?:.*\n)*?\nLA - .*?\n(?:.*\n)*?\nAID - .*?\n(?:.*\n)*?\n(?:AD  - (.*?@.*?\.\w+))', article)
            corresponding_email = email_match.group(1) if email_match else "Not Available"
            
            if company_affiliations:
                papers.append([pubmed_id, title, publication_date, ", ".join(non_academic_authors), ", ".join(company_affiliations), corresponding_email])
            
            if debug:
                print(f"Processed {pubmed_id}: {title}")
        
        return papers
    except Exception as e:
        print(f"Error fetching data: {e}", file=sys.stderr)
        return []

def save_to_csv(papers, filename):
    """Saves paper details to a CSV file."""
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["PubMed ID", "Title", "Publication Date", "Non-academic Authors", "Company Affiliations", "Corresponding Author Email"])
        writer.writerows(papers)
    print(f"Results saved to {filename}")

def main():
    parser = argparse.ArgumentParser(description="Fetch research papers from PubMed with company affiliations.")
    parser.add_argument("query", help="PubMed search query")  # ✅ Fixed: No duplicate argument
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug mode")
    parser.add_argument("-f", "--file", type=str, help="Output CSV file")

    args = parser.parse_args()
    
    papers = fetch_pubmed_papers(args.query, args.debug)
    
    if args.file:
        save_to_csv(papers, args.file)
    else:
        for paper in papers:
            print(" | ".join(paper))

if __name__ == "__main__":
    main()
