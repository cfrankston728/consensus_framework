"""
General Consensus Framework for Robust Statistical Locus Selection.

This package implements the mathematical framework described in:
"A General Consensus Framework for Robust Statistical Locus Selection"
by Frankston et al. (2025)

The framework integrates diverse sets of identified genomic loci (BED files)
into a robust consensus based on statistical dependencies among selectors.
"""

__version__ = "0.1.0"

# Import main functionality to expose at package level
from consensus_framework.locus_consensus import (
    load_basis_identifiers,
    create_prior_localizations,
    compute_identifier_selectivity,
    compute_consensus_statistics,
    simulate_independence_null_distribution,
    compute_p_value,
    apply_multiple_testing_correction,
    save_consensus_results,
)

__all__ = [
    "load_basis_identifiers",
    "create_prior_localizations",
    "compute_identifier_selectivity",
    "compute_consensus_statistics",
    "simulate_independence_null_distribution",
    "compute_p_value",
    "apply_multiple_testing_correction",
    "save_consensus_results",
]