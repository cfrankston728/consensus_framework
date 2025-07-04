#!/usr/bin/env python
"""
Core functions for the General Consensus Framework for Robust Statistical Locus Selection.

This module implements the mathematical framework described in:
"A General Consensus Framework for Robust Statistical Locus Selection"
by Frankston et al. (2025)

It provides functions for:
- Loading and processing genomic loci from BED files
- Creating prior localizations (connected components of L(∨X))
- Computing selectivity for basis identifiers
- Calculating consensus statistics and p-values
- Applying multiple testing corrections
- Saving results to BED files

The framework operates on three main assumptions:
1. Saturation - The union of all selectors covers the target loci
2. Independence of Incompetent Selection (IIS) - False loci are selected independently
3. Independence of Incompetent Omission (IIO) - True loci are omitted independently
"""

import numpy as np
import pandas as pd
from pathlib import Path
from statsmodels.stats.multitest import multipletests
import logging
from typing import List, Tuple, Dict, Union, Optional, Any

# Configure logging
logger = logging.getLogger(__name__)


def fast_indices_to_binary_vector(indices: List[int], vector_length: int) -> np.ndarray:
    """
    Convert a list of indices into a binary vector.
    
    This function efficiently creates a binary vector where positions
    specified by indices are set to 1, and all others are 0.
    
    Args:
        indices: List of integer indices to set to 1
        vector_length: Length of the resulting vector
        
    Returns:
        numpy.ndarray: Binary vector of specified length
    """
    binary_vector = np.zeros(vector_length, dtype=int)
    binary_vector[indices] = 1
    return binary_vector


def is_all_zeros_identifier_vector_vectorized(identifier_matrix: np.ndarray) -> np.ndarray:
    """
    Check which rows of an identifier matrix are all zeros.
    
    This function is used to identify positions where no identifiers
    are selecting, which mark the boundaries between prior localizations.
    
    Args:
        identifier_matrix: NumPy 2D array of identifier vectors
        
    Returns:
        numpy.ndarray: Boolean array indicating which rows are all zeros
    """
    return (identifier_matrix < 0.5).all(axis=1)


def load_basis_identifiers(location: str, sep: str, table_columns: Optional[str], 
                           location_column: Union[int, str] = 0) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """
    Load the basis identifiers table and all referenced BED files.
    
    This function loads the main table that references individual basis
    identifiers (BED files) and then loads all those BED files into a
    combined DataFrame. In the manuscript's notation, this creates the
    subset X ⊂ A of basis identifiers.
    
    Args:
        location: Path to the input BED file set table
        sep: Separator used in the input table
        table_columns: Comma-separated string of column labels
        location_column: Column index or name containing file paths
        
    Returns:
        tuple: (selected_loci_df, basis_identifiers_df, num_basis_identifiers)
            - selected_loci_df: DataFrame with all BED entries combined
            - basis_identifiers_df: DataFrame with the basis identifier table
            - num_basis_identifiers: Number of basis identifiers
    """
    # Process table columns
    if table_columns:
        names = table_columns.split(',')
    
    # Normalize separator
    if sep == '\\t':
        sep = '\t'
        
    logger.info("Reading basis identifier table...")
    basis_identifiers_df = pd.read_csv(location, sep=sep, header=None, names=names)
    logger.info("Basis identifier table loaded with %d rows", len(basis_identifiers_df))
    logger.debug("Basis identifier table head:\n%s", basis_identifiers_df.head())
    
    num_basis_identifiers = len(basis_identifiers_df)
    
    # Load all referenced BED files
    basis_identifier_df_list = []
    for idx, row in basis_identifiers_df.iterrows():
        x_location = row[location_column]
        x_sep = row['sep']
        x_table_columns = row['table_columns']
        
        logger.info("Loading loci data from: %s", x_location)
        
        if x_sep == '\\t':
            x_sep = '\t'
            
        try:
            df = pd.read_csv(
                x_location, 
                sep=x_sep, 
                header=None, 
                names=x_table_columns.split(','), 
                engine='python'
            )
            
            # Validate and clean BED data
            if 'start' in df.columns and 'end' in df.columns:
                # Remove invalid coordinates
                invalid_mask = (df['start'] < 0) | (df['end'] <= df['start'])
                if invalid_mask.any():
                    logger.warning(
                        "Removing %d invalid entries from %s (negative start or end <= start)",
                        invalid_mask.sum(), x_location
                    )
                    df = df[~invalid_mask]
            
            logger.debug("Loci data loaded for identifier_id %d with %d rows", idx, len(df))
            
            df['identifier_id'] = idx
            basis_identifier_df_list.append(df)
        except Exception as e:
            logger.error("Error loading %s: %s", x_location, str(e))
            raise
    
    # Combine all loaded BED files
    selected_loci_df = pd.concat(basis_identifier_df_list, ignore_index=True)
    logger.info("Combined loci DataFrame created with %d rows", len(selected_loci_df))
    logger.debug("Combined loci DataFrame head:\n%s", selected_loci_df.head())
    
    return selected_loci_df, basis_identifiers_df, num_basis_identifiers


def create_prior_localizations(selected_loci_df: pd.DataFrame, 
                              num_basis_identifiers: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create prior localizations from the combined loci DataFrame.
    
    This function implements the concept of saturation from the manuscript
    by creating a union of all regions selected by any identifier.
    
    The algorithm works by:
    1. Identifying all unique boundaries to create atomic regions
    2. Determining which identifiers are active in each atom
    3. Combining consecutive atoms where at least one identifier is active
       to form prior localizations (connected components of L(∨X))
    
    Args:
        selected_loci_df: DataFrame with all BED entries combined
        num_basis_identifiers: Number of basis identifiers
        
    Returns:
        tuple: (prior_localization_df, atoms_df)
            - prior_localization_df: DataFrame with prior localization regions in BED format
            - atoms_df: DataFrame with atomic regions and identifier vectors
    """
    logger.info("Creating prior localizations from %d loci", len(selected_loci_df))
    
    # Group by genomic coordinates and collect unique identifier IDs
    grouped_df = selected_loci_df.groupby(['chromosome', 'start', 'end'], as_index=False).agg({
        'identifier_id': lambda x: list(set(x))  # Collect all unique identifier IDs
    })
    logger.debug("Grouped DataFrame created with %d rows", len(grouped_df))
    
    # Convert identifier IDs to binary vectors
    grouped_df['differential_identifier_vector'] = grouped_df['identifier_id'].apply(
        lambda x: fast_indices_to_binary_vector(x, num_basis_identifiers)
    )
    grouped_df = grouped_df.drop(columns=['identifier_id'])
    
    # Create a copy of the dataframe and swap 'start' and 'end'
    # This is to handle interval boundaries properly
    swapped_df = grouped_df.copy()
    swapped_df[['start', 'end']] = swapped_df[['end', 'start']]
    swapped_df['differential_identifier_vector'] = swapped_df['differential_identifier_vector'].apply(lambda x: -1 * x)
    
    # Stack and sort the data
    stacked_df = pd.concat([grouped_df, swapped_df], ignore_index=True)
    atoms_df = stacked_df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
    atoms_df = atoms_df[['chromosome', 'start', 'differential_identifier_vector']]
    atoms_df = atoms_df.rename(columns={'start': 'position'})
    
    # Calculate true identifier vectors through integration (cumulative sum) of differentials
    # This creates atoms - regions where the identifier vector is constant
    atoms_df['identifier_vector'] = atoms_df['differential_identifier_vector'].cumsum()
    atoms_df['next_position'] = atoms_df['position'].shift(-1)
    # Handle the last position safely
    atoms_df.loc[atoms_df.index[-1], 'next_position'] = atoms_df.loc[atoms_df.index[-1], 'position'] + 1
    atoms_df['next_position'] = atoms_df['next_position'].astype(int)
    
    logger.debug("Created atoms DataFrame with %d rows", len(atoms_df))
    
    # Convert identifier vectors to a NumPy array for vectorized operations
    identifier_matrix = np.vstack(atoms_df['identifier_vector'].values)
    
    # Identify zero vector positions to determine prior localization boundaries
    # A zero vector means no identifiers are active, marking the boundary between prior localizations
    is_zero_vector = is_all_zeros_identifier_vector_vectorized(identifier_matrix)
    zero_indices = atoms_df.index[is_zero_vector]
    
    logger.debug("Found %d zero vector positions (boundaries between prior localizations)", len(zero_indices))
    
    # Construct prior localization regions by combining consecutive atoms with at least one identifier
    prior_localization_regions = []
    current_start = atoms_df.loc[atoms_df.index[0], 'position']
    
    for idx in zero_indices:
        prior_localization_regions.append([
            atoms_df.loc[idx, 'chromosome'], 
            current_start, 
            atoms_df.loc[idx, 'position']
        ])
        current_start = atoms_df.loc[idx, 'next_position']
    
    # Create DataFrame with prior localization regions
    prior_localization_df = pd.DataFrame(prior_localization_regions, columns=['chromosome', 'start', 'end'])
    logger.info("Created %d prior localization regions", len(prior_localization_df))
    logger.debug("Prior localization regions head:\n%s", prior_localization_df.head())
    
    return prior_localization_df, atoms_df


def compute_identifier_selectivity(atoms_df: pd.DataFrame, 
                                 num_basis_identifiers: int) -> Tuple[pd.DataFrame, np.ndarray, List, np.ndarray, int]:
    """
    Compute selectivity for each basis identifier.
    
    Selectivity is defined in the manuscript Section 2.3 as:
    px,K := EK[μ(L(x) ∩ K)/μ(K)]
    
    It is the expected probability that a given locus in a random prior 
    localization is identified by identifier x.
    
    Args:
        atoms_df: Sorted DataFrame with atomic regions and identifier vectors
        num_basis_identifiers: Number of basis identifiers
        
    Returns:
        tuple: (updated_atoms_df, identifier_coverage_sums, prior_localizations, selectivity_values, n_unique_instances)
            - updated_atoms_df: DataFrame with union identifier appended
            - identifier_coverage_sums: Sum of coverage for each identifier
            - prior_localizations: List of prior localizations
            - selectivity_values: Average coverage rates (selectivity)
            - n_unique_instances: Number of unique identifier vector patterns
    """
    logger.info("Computing identifier selectivity")
    
    # Extract identifier vectors as a matrix for vectorized operations
    identifier_matrix = np.vstack(atoms_df['identifier_vector'].values)
    
    # Create a new column with identifier vectors plus union
    union_indicator = (identifier_matrix > 0.5).any(axis=1).astype(int).reshape(-1, 1)
    identifier_matrix_with_union = np.hstack((identifier_matrix, union_indicator))
    
    # Add back to the DataFrame
    atoms_df['identifier_vector_with_union'] = [row for row in identifier_matrix_with_union]
    
    logger.debug("Added union identifier to atoms DataFrame")
    
    # Count unique identifier vector patterns (atomic patterns)
    # Convert to tuples for hashability in nunique
    identifier_vector_tuples = [tuple(row) for row in identifier_matrix_with_union]
    n_unique_instances = len(set(identifier_vector_tuples)) - 1  # Subtract 1 for all-zeros
    
    logger.info("Found %d unique atomic identifier-vectors", n_unique_instances)
    
    # Calculate identifier coverage sums across prior localizations
    vector_length = identifier_matrix_with_union.shape[1]
    identifier_coverage_sums = np.zeros(vector_length)
    current_prior_localization_start = atoms_df.loc[atoms_df.index[0], 'position']
    prior_localizations = []
    
    summed_vector = np.zeros(vector_length)
    
    # Process each row to compute coverage within prior localizations
    for idx, row in atoms_df.iterrows():
        current_identifier_vector = np.array(row['identifier_vector_with_union'])
        region_length = row['next_position'] - row['position']
        summed_vector += current_identifier_vector * region_length
        
        # If all identifiers are negative in this region (except possibly the union identifier)
        # This marks the boundary between prior localizations
        if all(current_identifier_vector[:-1] < 0.5):
            prior_localization_length = summed_vector[-1]
            if prior_localization_length > 0:
                coverage_vector = summed_vector / prior_localization_length
                identifier_coverage_sums += coverage_vector
            prior_localizations.append((current_prior_localization_start, row['next_position'], summed_vector.copy()))
            current_prior_localization_start = row['next_position']
            summed_vector = np.zeros(vector_length)
    
    # Handle any remaining region after the loop
    if current_prior_localization_start != atoms_df.loc[atoms_df.index[-1], 'next_position']:
        prior_localization_length = summed_vector[-1]
        if prior_localization_length > 0:
            coverage_vector = summed_vector / prior_localization_length
            identifier_coverage_sums += coverage_vector
        prior_localizations.append((current_prior_localization_start, atoms_df.loc[atoms_df.index[-1], 'next_position'], summed_vector.copy()))
    
    # Calculate average coverage (selectivity) across prior localizations
    selectivity_values = identifier_coverage_sums / len(prior_localizations)
    logger.info("Computed selectivity values for %d identifiers", len(selectivity_values) - 1)
    logger.debug("Selectivity values: %s", selectivity_values[:-1])  # Exclude union
    
    return atoms_df, identifier_coverage_sums, prior_localizations, selectivity_values, n_unique_instances


def compute_consensus_statistics(identifier_vector_with_union: np.ndarray, 
                                selectivity_values: np.ndarray) -> Tuple[float, float]:
    """
    Compute consensus selection and omission statistics for a given identifier vector.
    
    This implements the key information-theoretic statistics from Section 2.6 of the manuscript:
    - Consensus selection statistic: 1/|X| * sum(log(1/px,K)) for selected loci
    - Consensus omission statistic: 1/|X| * sum(log(1/(1-px,K))) for omitted loci
    
    Args:
        identifier_vector_with_union: Binary vector indicating which identifiers selected this locus,
                                    including the union identifier at the end
        selectivity_values: Vector of selectivity values for each identifier
        
    Returns:
        tuple: (consensus_selection_score, consensus_omission_score)
    """
    # Remove the union identifier for these calculations
    identifier_vector = identifier_vector_with_union[:-1]
    selectivities = selectivity_values[:-1]  # Exclude the union identifier
    
    # Identify selection and omission indices using boolean masks
    selection_mask = identifier_vector > 0.5
    omission_mask = ~selection_mask
    
    # Calculate consensus scores
    # Use masked arrays for cleaner calculation
    consensus_selection_score = np.sum(-np.log(selectivities[selection_mask]))
    consensus_omission_score = np.sum(-np.log(1 - selectivities[omission_mask]))
    
    return consensus_selection_score, consensus_omission_score


def simulate_independence_null_distribution(selectivity_values: np.ndarray, 
                                          n_simulations: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate null distributions for consensus selection and omission statistics.
    
    This creates empirical null distributions by simulating random identifier
    selections based on the observed selectivity rates, assuming independence.
    This directly tests the independence of incompetent selection (IIS) and 
    independence of incompetent omission (IIO) assumptions from Section 2.4.
    
    This implementation ensures exactly n_simulations samples where at least
    one identifier selects, using adaptive oversampling.
    
    Args:
        selectivity_values: Vector of selectivity values for each identifier
        n_simulations: Number of simulations to perform
        
    Returns:
        tuple: (sorted_selection_stats, selection_ecdf, sorted_omission_stats, omission_ecdf)
            - sorted_selection_stats: Sorted selection statistics from simulation
            - selection_ecdf: Empirical CDF for selection statistics
            - sorted_omission_stats: Sorted omission statistics from simulation
            - omission_ecdf: Empirical CDF for omission statistics
    """
    logger.info("Simulating null distribution with %d simulations", n_simulations)
    
    # Use only non-union selectivities
    selectivities = selectivity_values[:-1]
    R = len(selectivities)
    
    # Calculate probability of all zeros (no identifiers selecting)
    p_none = np.prod(1 - selectivities)
    logger.debug("Probability of no identifiers selecting: %.4f", p_none)
    
    # If this probability is high, we'll need substantial oversampling
    if p_none > 0.1:  # Threshold can be adjusted
        # Initial estimate of samples needed with 20% safety margin
        oversample_factor = 1 / (1 - p_none) * 1.2
        target_samples = int(n_simulations * oversample_factor)
        logger.debug("High zero probability, using oversample factor: %.2f", oversample_factor)
        
        # Start with this estimate, then adapt if needed
        sims = np.random.binomial(1, p=selectivities, size=(target_samples, R))
        non_zero_sims = sims[sims.sum(axis=1) > 0]
        
        # If we didn't get enough samples, incrementally generate more
        while len(non_zero_sims) < n_simulations:
            logger.debug("Insufficient samples (%d/%d), generating more", len(non_zero_sims), n_simulations)
            # Add another batch with the same size as the first
            additional = np.random.binomial(1, p=selectivities, size=(target_samples, R))
            additional_non_zero = additional[additional.sum(axis=1) > 0]
            non_zero_sims = np.vstack((non_zero_sims, additional_non_zero))
        
        # Take exactly n_simulations samples
        sims = non_zero_sims[:n_simulations]
    else:
        # If probability of zeros is low, modest oversampling is sufficient
        target_samples = int(n_simulations * 1.2)  # 20% oversampling
        logger.debug("Low zero probability, using standard 20%% oversampling")
        sims = np.random.binomial(1, p=selectivities, size=(target_samples, R))
        non_zero_sims = sims[sims.sum(axis=1) > 0]
        
        # Ensure we have exactly n_simulations
        if len(non_zero_sims) >= n_simulations:
            sims = non_zero_sims[:n_simulations]
        else:
            # Unlikely case, but handle it
            logger.debug("Insufficient samples, adding more in batches")
            while len(non_zero_sims) < n_simulations:
                additional = np.random.binomial(1, p=selectivities, size=(n_simulations, R))
                additional_non_zero = additional[additional.sum(axis=1) > 0]
                non_zero_sims = np.vstack((non_zero_sims, additional_non_zero))
            sims = non_zero_sims[:n_simulations]
    
    logger.debug("Generated %d simulation samples", len(sims))
    
    # Compute stats
    sel_weights = -np.log(selectivities)
    omit_weights = -np.log(1 - selectivities)
    
    selection_stats = sims.dot(sel_weights)
    omission_stats = (1 - sims).dot(omit_weights)
    
    # Sort & build ECDF
    sorted_selection_stats = np.sort(selection_stats)
    sorted_omission_stats = np.sort(omission_stats)
    selection_ecdf = np.arange(1, n_simulations+1) / n_simulations
    omission_ecdf = selection_ecdf.copy()
    
    return sorted_selection_stats, selection_ecdf, sorted_omission_stats, omission_ecdf


def compute_p_value(obs: float, sorted_stats: np.ndarray, ecdf: np.ndarray) -> float:
    """
    Compute p-value using the empirical cumulative distribution function.
    
    This calculates 1 - ECDF(observed value), representing the probability 
    of observing a statistic as or more extreme than the observed value
    under the null hypothesis of independence.
    
    Args:
        obs: Observed statistic value
        sorted_stats: Sorted statistics from null distribution
        ecdf: Empirical CDF corresponding to sorted_stats
        
    Returns:
        float: P-value (1 - ECDF at the observed value)
    """
    # Find the index where the observed value would be inserted
    idx = np.searchsorted(sorted_stats, obs, side='right')
    return 1 - ecdf[idx-1] if idx > 0 else 1.0


def apply_multiple_testing_correction(atoms_df: pd.DataFrame, 
                                    n_unique_instances: int, 
                                    alpha: float, 
                                    mtc: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Apply multiple testing correction to the p-values.
    
    This function implements various multiple testing correction methods
    using the statsmodels.stats.multitest.multipletests function.
    
    Args:
        atoms_df: DataFrame with p-values
        n_unique_instances: Number of hypothesis tests (m)
        alpha: Significance threshold
        mtc: Multiple testing correction method ('raw', 'bonf', 'sidak', 'bh')
        
    Returns:
        tuple: (corrected_df, consensus_regions_sel, consensus_regions_omit)
            - corrected_df: DataFrame with corrected p-values
            - consensus_regions_sel: Significant regions based on selection
            - consensus_regions_omit: Significant regions based on omission
    """
    logger.info("Applying %s multiple testing correction with alpha=%.4f", mtc, alpha)
    
    # Map mtc to statsmodels method names
    mtc_method_map = {
        'raw': None,  # No correction
        'bonf': 'bonferroni',
        'sidak': 'sidak',
        'bh': 'fdr_bh'
    }
    
    # Number of hypothesis tests
    m = n_unique_instances
    logger.debug("Number of tests (m): %d", m)
    
    consensus_regions = {}
    
    # Apply corrections for both selection and omission statistics
    for stat, stat_name in [('selection', 'sel'), ('omission', 'omit')]:
        logger.debug("Processing %s statistics", stat)
        
        # Extract p-values for this statistic
        p_values = atoms_df[f'p_{stat}'].values
        
        # Apply correction if method is specified
        if mtc != 'raw':
            method = mtc_method_map[mtc.lower()]
            _, corrected_pvals, _, _ = multipletests(
                p_values, 
                alpha=alpha, 
                method=method
            )
            atoms_df[f'p_{stat}_{mtc.lower()}'] = corrected_pvals
            selection_col = f'p_{stat}_{mtc.lower()}'
        else:
            selection_col = f'p_{stat}'
        
        # Flag significant bins
        atoms_df['is_significant'] = atoms_df[selection_col] <= alpha
        atoms_df['block_change'] = atoms_df['is_significant'] != atoms_df['is_significant'].shift(1, fill_value=False)
        atoms_df['block_id'] = atoms_df['block_change'].cumsum()
        
        # Group & collapse into regions
        sig_blocks = atoms_df[atoms_df['is_significant']].groupby('block_id')
        
        if len(sig_blocks) > 0:
            bed_regions = sig_blocks.agg({
                'chromosome': 'first',
                'position': 'first',
                'next_position': 'last'
            }).reset_index(drop=True)
            
            # Rename columns to match BED format
            bed_regions = bed_regions.rename(columns={'position': 'start', 'next_position': 'end'})
            
            logger.info("Found %d significant %s regions", len(bed_regions), stat)
            logger.debug("Significant regions head:\n%s", bed_regions.head() if len(bed_regions) > 0 else "None")
        else:
            bed_regions = pd.DataFrame(columns=['chromosome', 'start', 'end'])
            logger.info("No significant %s regions found", stat)
        
        # Store the results
        consensus_regions[stat_name] = bed_regions
    
    return atoms_df, consensus_regions['sel'], consensus_regions['omit']


def save_consensus_results(output_dir: Union[str, Path], 
                           prior_localization_regions_df: pd.DataFrame, 
                           consensus_track_df: pd.DataFrame, 
                           consensus_regions_sel: Optional[pd.DataFrame] = None,
                           consensus_regions_omit: Optional[pd.DataFrame] = None,
                           alpha: Optional[float] = None, 
                           mtc: Optional[str] = None) -> None:
    """
    Save the consensus results to output files.
    
    Args:
        output_dir: Directory to save output files
        prior_localization_regions_df: DataFrame with prior localization regions
        consensus_track_df: DataFrame with consensus track
        consensus_regions_sel: DataFrame with significant regions based on selection
        consensus_regions_omit: DataFrame with significant regions based on omission
        alpha: Significance threshold
        mtc: Multiple testing correction method
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    logger.info("Saving results to %s", output_dir)
    
    # Save prior localization regions
    prior_localization_regions_file = output_dir / 'prior_localization_regions.bed'
    prior_localization_regions_df.to_csv(prior_localization_regions_file, sep='\t', header=False, index=False)
    logger.info("Prior localization regions saved to %s", prior_localization_regions_file)
    
    # Save consensus track
    consensus_track_file = output_dir / 'consensus_track.bed'
    consensus_track_df.to_csv(consensus_track_file, sep='\t', header=True, index=False)
    logger.info("Consensus track with p-values saved to %s", consensus_track_file)
    
    # Save significant regions if available
    if alpha is not None:
        if consensus_regions_sel is not None and not consensus_regions_sel.empty:
            sel_file = output_dir / f'consensus_regions.selection.alpha-{alpha}.mtc-{mtc}.bed'
            consensus_regions_sel.to_csv(sel_file, sep='\t', header=False, index=False)
            logger.info("Selection regions saved to %s", sel_file)
            
        if consensus_regions_omit is not None and not consensus_regions_omit.empty:
            omit_file = output_dir / f'consensus_regions.omission.alpha-{alpha}.mtc-{mtc}.bed'
            consensus_regions_omit.to_csv(omit_file, sep='\t', header=False, index=False)
            logger.info("Omission regions saved to %s", omit_file)