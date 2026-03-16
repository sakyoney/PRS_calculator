# PRS Calculator Web App

A minimal Flask-based PRS calculator with a simple HTML interface.

## Features

- Upload GWAS and sample files from the browser
- Supports CSV, TSV, TXT, XLS, XLSX
- Accepts either genotype or dosage in the sample file
- Displays PRS score in the browser
- Lets you download the full SNP-level result as CSV

## Expected input

### GWAS file
Required columns (case-insensitive matching):
- `SNP`
- `BETA`
- `EFFECT_ALLELE`

### Sample file
Required:
- `SNP`

And one of:
- `DOSAGE`
- `GENOTYPE`

Genotype examples:
- `A/G`
- `AA`
- `C|T`

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Then open:

```text
http://127.0.0.1:5000
```

## Example files

- `example_gwas.tsv`
- `example_sample_genotype.csv`
- `example_sample_dosage.csv`

## Notes

This is a lightweight demo/MVP. For production or research-grade use, you would likely want to add:
- allele harmonization
- strand handling
- sample identifiers / multi-sample support
- validation and logging
- async processing for large files
