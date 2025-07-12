import re
from typing import Union

# Keywords used to detect pharmaceutical/biotech companies
COMPANY_KEYWORDS = [
    "pharma", "biotech", "therapeutics", "inc", "ltd",
    "gmbh", "corp", "company", "llc", "diagnostic"
]

def is_non_academic(affiliation: str) -> bool:
    """
    Heuristic to detect non-academic institutions (i.e., companies).
    Filters out universities, colleges, hospitals, etc.
    """
    affiliation = affiliation.lower()
    return (
        any(keyword in affiliation for keyword in COMPANY_KEYWORDS)
        and not any(acad in affiliation for acad in ["university", "college", "institute", "hospital", "lab"])
    )

def extract_email(text: Union[str, list]) -> str:
    """Extract the first valid email from a text or list of strings."""
    if isinstance(text, list):
        text = " ".join(text)
    match = re.search(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", text)
    return match.group(0) if match else ""
