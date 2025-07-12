# 🧬 PubMed Paper Fetcher CLI

A Python CLI tool to fetch research papers from PubMed based on a user-defined query, and filter results with at least one non-academic author affiliated with a pharmaceutical or biotech company.

---

## 📌 Features

- 🔍 Query PubMed using full-text search
- 🧑‍🔬 Detect non-academic authors using email/affiliation heuristics
- 🏢 Identify authors affiliated with pharma/biotech companies
- 📤 Output results to console or CSV
- ⚙️ Supports CLI flags: `--debug`, `--file`
- 🧪 Built with Poetry and Biopython

---

## 🚀 Installation

### 1. Clone the repository

```bash
git clone https://github.com/Amal047/Python-program-to-fetch-research-papers-based-on-a-user-specified-query.git
cd Python-program-to-fetch-research-papers-based-on-a-user-specified-query
```

### 2. Install using Poetry

Make sure [Poetry](https://python-poetry.org/docs/#installation) is installed:

```bash
poetry install
```

---

## 🛠️ Usage

### Basic usage

```bash
poetry run get-papers-list "diabetes treatment"
```

### Save results to a CSV

```bash
poetry run get-papers-list "cancer immunotherapy" --file results.csv
```

### Enable debug output

```bash
poetry run get-papers-list "AI in pharma" --debug
```

---

## 📤 Output Format

CSV and console output includes:

- **PubmedID**
- **Title**
- **Publication Date**
- **Non-academic Author(s)**
- **Company Affiliation(s)**
- **Corresponding Author Email**

---

## 🧠 Heuristics for Company Detection

An author is considered **non-academic** if their affiliation:
- Contains keywords like `pharma`, `biotech`, `therapeutics`, `inc`, `corp`, `ltd`, etc.
- **Does not** include words like `university`, `institute`, `college`, `hospital`, or `lab`

Emails are extracted using regex from affiliation strings.

---

## 🧱 Project Structure

```
get_papers/
├── cli.py         # CLI entry point
├── fetch.py       # PubMed API functions
├── filter.py      # Filtering logic for companies/emails
├── output.py      # Console/CSV output functions
```

---

## 🧪 Dependencies

- [Poetry](https://python-poetry.org/)
- [Biopython](https://biopython.org/)
- [Requests](https://docs.python-requests.org/)
- [Tqdm](https://github.com/tqdm/tqdm)

Install with:
```bash
poetry install
```

---

## ✅ Tasks Completed (per assignment)

- [x] PubMed API integration
- [x] Filtering based on non-academic affiliations
- [x] CSV and console output
- [x] CLI with flags `--file`, `--debug`
- [x] Poetry-based project setup
- [x] Git version control on GitHub
- [x] Typed Python (using type hints)
- [x] Modular code structure
- [x] README with usage instructions
- [ ] (Optional) Publish on TestPyPI

---


## 📦 TestPyPI Package

## This package is published to [TestPyPI](https://test.pypi.org/) for demonstration purposes:

```bash
 pip install --index-url https://test.pypi.org/simple/ get-papers-nikhilvh
```


## 👨‍💻 Author

**Nikhil V**  
📧 nikhilvh2600@gmail.com

---

## 📄 License

This project is licensed under the MIT License.
