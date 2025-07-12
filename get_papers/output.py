import csv
from typing import List, Dict

def write_to_csv(filename: str, papers: List[Dict]):
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["PubmedID", "Title", "Publication Date", "Non-academic Author(s)", "Company Affiliation(s)", "Corresponding Author Email"])
        writer.writeheader()
        writer.writerows(papers)

def print_to_console(papers: List[Dict]):
    for paper in papers:
        print("\n---")
        for key, value in paper.items():
            print(f"{key}: {value}")
