# Consensus Framework for Robust Statistical Locus Selection

This package implements the mathematical framework described in "A General Consensus Framework for Robust Statistical Locus Selection" by Connor M. Frankston.

## Overview

The consensus framework integrates diverse sets of genomic loci (from multiple BED files) into a robust consensus based on statistical dependencies among selectors. It operates on three main assumptions:

1. **Saturation** - The union of all selectors covers the target loci
2. **Independence of Incompetent Selection (IIS)** - False loci are selected independently
3. **Independence of Incompetent Omission (IIO)** - True loci are omitted independently

## Installation

```bash
# Clone the repository
git clone https://github.com/frankston/consensus-framework.git
cd consensus-framework

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

## Usage

### Command Line Interface

```bash
# Basic usage
consensus-framework --location input_table.tsv --output-dir ./results

# With multiple testing correction and significance threshold
consensus-framework --location input_table.tsv --output-dir ./results --alpha 0.05 --mtc bh

# Enable verbose logging
consensus-framework --location input_table.tsv --output-dir ./results --verbose
```

### Options

- `--location`: Path to the input table file containing BED file references (required)
- `--sep`: Separator used in the input table (default: tab)
- `--table-columns`: Comma-separated string of column labels
- `--output-dir`: Directory to save output files (required)
- `--skip-p-vals/--get-p-vals`: Skip or compute p-values (default: compute)
- `--alpha`: Significance threshold for consensus computation
- `--n-simulations`: Number of simulations for null distribution (default: 10000)
- `--recompute`: Recompute outputs even if existing files are found
- `--mtc`: Multiple-testing correction method (`raw`, `bonf`, `sidak`, `bh`) (default: `bh`)
- `--verbose`: Enable verbose (debug) logging

### Programmatic Usage

```python
from consensus_framework import (
    load_basis_identifiers,
    create_prior_localizations,
    compute_identifier_selectivity,
    apply_multiple_testing_correction
)

# Load identifiers
selected_loci_df, basis_identifiers_df, num_basis_identifiers = load_basis_identifiers(
    "input_table.tsv", "\t", "location,sep,table_columns"
)

# Create prior localizations
prior_localization_df, atoms_df = create_prior_localizations(
    selected_loci_df, num_basis_identifiers
)

# Compute selectivity
atoms_df, coverage_sums, prior_localizations, selectivity_values, n_unique = compute_identifier_selectivity(
    atoms_df, num_basis_identifiers
)

# ... additional processing steps
```

### Input Format

The input table should contain information about each BED file to be used as a basis identifier:

```
location    sep    table_columns
file1.bed    \t    chromosome,start,end
file2.bed    \t    chromosome,start,end
```

### Output Files

- `prior_localization_regions.bed`: BED file with prior localization regions (union of all identifiers)
- `consensus_track.bed`: BED file with consensus statistics and p-values
- `consensus_regions.selection.alpha-X.mtc-Y.bed`: Significant regions based on selection (if alpha is provided)
- `consensus_regions.omission.alpha-X.mtc-Y.bed`: Significant regions based on omission (if alpha is provided)

## Implementation Details

The framework is implemented in two main modules:

1. `locus_consensus.py`: Contains the core statistical consensus functionality
2. `locus_consensus_cli.py`: Provides the command-line interface

### Key Components

1. **Data Loading**: Load basis identifiers and their corresponding BED files
2. **Prior Localization Construction**: Create atomic regions and prior localizations
3. **Selectivity Computation**: Calculate selectivity for each basis identifier
4. **Consensus Statistics**: Compute selection and omission statistics
5. **Null Distribution Simulation**: Generate empirical null distributions
6. **Multiple Testing Correction**: Apply corrections (Bonferroni, Šidák, Benjamini-Hochberg)
7. **Result Saving**: Output consensus regions as BED files

## Mathematical Framework

The package implements the mathematical framework described in the manuscript, including:

- Creating prior localizations as connected components of the union of all identifiers
- Computing selectivity as the expected probability of identifier selection
- Calculating consensus statistics based on information-theoretic principles
- Testing against independence hypotheses through simulation
- Applying multiple testing correction to control error rates

## Testing

Run the unit tests with pytest:

```bash
pytest tests/
```

## Performance Considerations

- The implementation uses vectorized operations where possible for improved performance
- For large datasets, consider adjusting the `--n-simulations` parameter to balance accuracy and speed
- Memory usage scales with the number and size of input BED files

## License

[MIT License](LICENSE)
