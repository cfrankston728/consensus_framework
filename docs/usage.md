# Consensus Framework Usage Guide

This document provides detailed usage instructions for the Consensus Framework for Robust Statistical Locus Selection.

## Command Line Interface

The consensus framework provides a command-line interface for processing BED files and generating consensus regions.

### Basic Usage

```bash
consensus-framework --location input_table.tsv --output-dir ./results
```

### Required Arguments

- `--location`: Path to the input table file containing BED file references
- `--output-dir`: Directory to save output files

### Optional Arguments

- `--sep`: Separator used in the input table (default: tab)
- `--table-columns`: Comma-separated string of column labels
- `--skip-p-vals/--get-p-vals`: Skip or compute p-values (default: compute)
- `--alpha`: Significance threshold for consensus computation
- `--n-simulations`: Number of simulations for null distribution (default: 10000)
- `--recompute`: Recompute outputs even if existing files are found
- `--mtc`: Multiple-testing correction method (default: `bh`)
  - `raw`: Uncorrected p-values
  - `bonf`: Bonferroni correction
  - `sidak`: Šidák correction
  - `bh`: Benjamini-Hochberg (FDR) correction
- `--verbose`: Enable verbose (debug) logging

### Examples

#### Basic Consensus Generation

```bash
consensus-framework --location input_table.tsv --output-dir ./results
```

#### With Significance Threshold and Multiple Testing Correction

```bash
consensus-framework --location input_table.tsv --output-dir ./results --alpha 0.05 --mtc bh
```

#### Skip P-value Computation (Faster)

```bash
consensus-framework --location input_table.tsv --output-dir ./results --skip-p-vals
```

#### With Verbose Logging

```bash
consensus-framework --location input_table.tsv --output-dir ./results --verbose
```

#### With Custom Input Format

```bash
consensus-framework --location input_table.csv --sep "," --table-columns "file_path,separator,column_names" --output-dir ./results
```

## Input Format

The input table should contain information about each BED file to be used as a basis identifier. The default format is a tab-separated file with the following columns:

```
location    sep    table_columns
file1.bed    \t    chromosome,start,end
file2.bed    \t    chromosome,start,end
```

- `location`: Path to the BED file
- `sep`: Separator used in the BED file (use `\t` for tab)
- `table_columns`: Comma-separated list of column names in the BED file

## Output Files

The framework generates several output files:

1. `prior_localization_regions.bed`: BED file with prior localization regions (union of all identifiers)
2. `consensus_track.bed`: BED file with consensus statistics and p-values, with columns:
   - `chromosome`: Chromosome name
   - `start`: Start position
   - `end`: End position
   - `consensus_selection_score`: Information-theoretic score for selection
   - `p_selection`: P-value for the selection score
   - `consensus_omission_score`: Information-theoretic score for omission
   - `p_omission`: P-value for the omission score

3. If a significance threshold (`--alpha`) is provided:
   - `consensus_regions.selection.alpha-X.mtc-Y.bed`: Significant regions based on selection
   - `consensus_regions.omission.alpha-X.mtc-Y.bed`: Significant regions based on omission

## Programmatic Usage

You can also use the consensus framework programmatically in your Python code:

```python
import logging
from pathlib import Path
from consensus_framework import (
    load_basis_identifiers,
    create_prior_localizations,
    compute_identifier_selectivity,
    compute_consensus_statistics,
    simulate_independence_null_distribution,
    apply_multiple_testing_correction,
    save_consensus_results
)

# Configure logging
logging.basicConfig(level=logging.INFO)

# Load basis identifiers
selected_loci_df, basis_identifiers_df, num_basis_identifiers = load_basis_identifiers(
    "input_table.tsv", "\t", "location,sep,table_columns"
)

# Create prior localizations
prior_localization_df, atoms_df = create_prior_localizations(
    selected_loci_df, num_basis_identifiers
)

# Compute selectivity
atoms_df, identifier_coverage_sums, prior_localizations, selectivity_values, n_unique_instances = compute_identifier_selectivity(
    atoms_df, num_basis_identifiers
)

# Update basis identifiers with selectivity
basis_identifiers_df['selectivity'] = selectivity_values[:-1]  # Exclude union

# Ensure no zero selectivities
selectivities = np.clip(selectivity_values, 1e-10, None)

# Extract identifier vectors and compute consensus statistics
# [Additional code to compute consensus statistics]

# Simulate null distributions
sorted_selection_stats, selection_ecdf, sorted_omission_stats, omission_ecdf = simulate_independence_null_distribution(
    selectivities, 10000
)

# Compute p-values
# [Additional code to compute p-values]

# Apply multiple testing correction
atoms_df, consensus_regions_sel, consensus_regions_omit = apply_multiple_testing_correction(
    atoms_df, n_unique_instances, 0.05, 'bh'
)

# Save results
output_dir = Path("./results")
save_consensus_results(
    output_dir=output_dir,
    prior_localization_regions_df=prior_localization_df,
    consensus_track_df=track_df,
    consensus_regions_sel=consensus_regions_sel,
    consensus_regions_omit=consensus_regions_omit,
    alpha=0.05,
    mtc='bh'
)
```

## Performance Tips

- For large datasets, reduce the number of simulations (e.g., `--n-simulations 1000`) for faster processing
- Use the `--skip-p-vals` option if you only need the prior localization regions
- Consider filtering input BED files to focus on specific regions of interest
- The memory usage scales with the number and size of input BED files