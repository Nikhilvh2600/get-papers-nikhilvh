import argparse
from typing import List, Dict
from get_papers.fetch import search_pubmed, fetch_details
from get_papers.filter import is_non_academic, extract_email
from get_papers.output import write_to_csv, print_to_console
from Bio import Medline


def parse_articles(raw_data: str, debug: bool = False) -> List[Dict[str, str]]:
    """Parse MEDLINE-formatted raw data and extract relevant paper info."""
    records = list(Medline.parse(raw_data.splitlines()))
    results = []

    for record in records:
        pmid = record.get("PMID", "")
        title = record.get("TI", "")
        date = record.get("DP", "")
        authors = record.get("AU", [])
        affiliations = record.get("AD", "")

        # Handle cases where affiliations are a list
        if isinstance(affiliations, list):
            affiliations = " | ".join(affiliations)

        email = extract_email(affiliations)

        if is_non_academic(affiliations):
            results.append({
                "PubmedID": pmid,
                "Title": title,
                "Publication Date": date,
                "Non-academic Author(s)": "; ".join(authors),
                "Company Affiliation(s)": affiliations,
                "Corresponding Author Email": email
            })
        elif debug:
            print(f"Skipping academic paper: {title}")

    return results


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Fetch PubMed papers with non-academic authors."
    )
    parser.add_argument("query", type=str, help="PubMed search query")
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Enable debug output"
    )
    parser.add_argument(
        "-f", "--file", type=str, help="Save results to CSV file"
    )

    args = parser.parse_args()
    if args.debug:
        print(f"Searching PubMed for: {args.query}")

    ids = search_pubmed(args.query)
    if args.debug:
        print(f"Found {len(ids)} paper IDs")

    raw_data = fetch_details(ids)
    parsed_data = parse_articles(raw_data, debug=args.debug)

    if args.file:
        write_to_csv(args.file, parsed_data)
        if args.debug:
            print(f"Results saved to {args.file}")
    else:
        print_to_console(parsed_data)


if __name__ == "__main__":
    main()
