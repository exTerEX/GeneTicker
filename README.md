
# GeneTicker

A small utility for calculating codon frequencies from GenBank (GBFF) files.

This project parses CDS features from GenBank flat files, counts codon occurrences,
maps codons to amino acids (including marking start/stop codons), and produces
tabular frequency reports that can be exported to several formats or printed to
the console.

## Features
- Parse GenBank (`.gbff`) files using Biopython
- Count codons across CDS features and compute frequency (percentage)
- Mark special codons as `START` or `STOP` when applicable
- Export results to `.csv`, `.tsv`, `.json`, `.pickle`, `.xlsx`, `.xml`,
	and columnar formats (feather/parquet/orc) if dependencies are available
- Simple CLI entry points and programmatic API

## Quickstart

1. Install GeneTicker:

Install a specific version by changing v0.1.0 to a specific tag version found on github.

```bash
pip install git+https://github.com/exTerEX/GeneTicker@v0.1.0
```

If you want to install optional dependencies (all, excel, hadoop, arrow, parquet, orc)

```bash
pip install GeneTicker[all]@git+https://github.com/exTerEX/GeneTicker
```

or

```bash
pip install GeneTicker[all]@git+https://github.com/exTerEX/GeneTicker@v0.1.0
```

2. Example usage:

**Command-line interface (CLI):**

```bash
gtick -o output.csv path/to/input.gbff
```

**Callable:**

```python
from pathlib import Path
from GeneTicker.core import run_codon_analysis

df = run_codon_analysis(Path('path/to/input.gbff'))
```

Output a Pandas DataFrame with columns:
- *Amino Acid*: Amino acid IUPAC one letter codes
- *Codon*: Nucleic acid three letter codes
- *Count*: Number of specific codon
- *Freq. (%)*: Frequency of total codon count
- *Special Type*: Indicate special role (START, STOP, etc.)

### Notes and dependencies

- Biopython is required to parse GenBank files (`Bio` package).
- Optional packages like `pyarrow`, `fastparquet`, or `lxml` are required to export to some formats (parquet/feather/orc/xml).

## License

This repository is provided as-is. See `pyproject.toml` for package metadata.

