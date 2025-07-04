#!/usr/bin/env python
"""
Command-line interface for the General Consensus Framework for Robust Statistical Locus Selection.

This module provides a command-line interface to the core functions in the consensus framework.
It handles argument parsing, logging configuration, and orchestration of the consensus computation.

Usage:
    python cli.py --location input_table.tsv --output-dir ./results
    python cli.py --location input_table.tsv --output-dir ./results --alpha 0.05 --mtc bh --verbose
"""

import click
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import sys
from typing import Optional

# Import core functions
from locus_consensus import (
    load_basis_identifiers,
    create_prior_localizations,
    compute_identifier_selectivity,
    compute_consensus_statistics,
    simulate_independence_null_distribution,
    compute_p_value,
    apply_multiple_testing_correction,
    save_consensus_results
)

# Configure root logger
logger = logging.getLogger()


def setup_logging(verbose: bool) -> None:
    """Configure logging based on verbosity level."""
    # Clear any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Set up console handler
    console_handler = logging.StreamHandler(sys.stdout)
    
    # Define format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    
    # Set log level based on verbosity
    if verbose:
        logger.setLevel(logging.DEBUG)
        console_handler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        console_handler.setLevel(logging.INFO)
    
    # Add handler to logger
    logger.addHandler(console_handler)


@click.command()
@click.option('--location', type=click.Path(exists=True), required=True, 
              help="Location of the input bed file set table file.")
@click.option('--sep', default='\t', 
              help="Separator for the input table file.")
@click.option('--table-columns', default=None, 
              help="String of column labels joined by commas.")
@click.option('--output-dir', type=click.Path(), required=True, 
              help="Directory to save output files.")
@click.option('--skip-p-vals/--get-p-vals', default=False, 
              help='Skip or compute p-values. Default is to compute p-values.')
@click.option('--alpha', default=None, type=float, 
              help='Significance threshold for consensus computation. If None, outputs the consensus track without thresholding.')
@click.option('--n-simulations', default=10000, 
              help='Number of simulations for null distribution.')
@click.option('--recompute', is_flag=True, default=False, 
              help='If set, recomputes outputs and overwrites existing files.')
@click.option('--mtc', type=click.Choice(['raw','bonf','sidak','bh'], case_sensitive=False), default='bh',
              help=("Which multiple-testing correction to include in the consensus track: "
                   "`raw` = uncorrected p_value; "
                   "`bonf` = Bonferroni; "
                   "`sidak` = Šidák; "
                   "`bh` = Benjamini–Hochberg (FDR)."))
@click.option('--verbose', is_flag=True, default=False,
              help="Enable verbose (debug) logging.")
def construct_consensus_beds(location: str, 
                            sep: str, 
                            table_columns: Optional[str], 
                            output_dir: str, 
                            skip_p_vals: bool, 
                            alpha: Optional[float], 
                            n_simulations: int, 
                            recompute: bool, 
                            mtc: str,
                            verbose: bool) -> None:
    """
    Construct consensus BED files from multiple basis identifiers.
    
    This command-line tool implements the consensus framework described in
    "A General Consensus Framework for Robust Statistical Locus Selection"
    by Frankston et al. (2025).
    
    The tool takes a table of BED files as input (basis identifiers), computes prior
    localizations (union of all regions), calculates statistical scores based on
    identifier selectivity, and identifies significant consensus regions using
    various statistical tests.
    
    Examples:
        python cli.py --location input_table.tsv --output-dir ./results
        python cli.py --location input_table.tsv --output-dir ./results --alpha 0.05 --mtc bh --verbose
    """
    # Configure logging
    setup_logging(verbose)
    
    # Convert output_dir to Path
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Define output file paths
    prior_localization_regions_file = output_path / 'prior_localization_regions.bed'
    consensus_track_file = output_path / 'consensus_track.bed'
    
    if alpha is not None:
        consensus_regions_file = output_path / f'consensus_regions.alpha-{alpha}.mtc-{mtc.lower()}.bed'
    
    # Initialize flags to determine whether to recompute files
    recompute_prior_localizations = recompute
    recompute_consensus_track = recompute
    
    get_p_vals = not skip_p_vals
    
    # Check if prior_localization_regions.bed exists and is valid
    if prior_localization_regions_file.exists() and not recompute:
        try:
            prior_localization_df = pd.read_csv(prior_localization_regions_file, sep='\t', header=None)
            logger.info("Using existing prior localization regions file: %s", prior_localization_regions_file)
        except Exception as e:
            logger.error("Error reading %s: %s", prior_localization_regions_file, str(e))
            logger.info("Backing up and recreating prior localization regions file")
            # Create backup before removing
            backup_file = prior_localization_regions_file.with_suffix('.bed.bak')
            prior_localization_regions_file.rename(backup_file)
            recompute_prior_localizations = True
    else:
        recompute_prior_localizations = True
    
    # Check if consensus_track.bed exists and is valid
    if get_p_vals and consensus_track_file.exists() and not recompute:
        try:
            track_df = pd.read_csv(consensus_track_file, sep='\t', header=0)
            logger.info("Using existing consensus track file: %s", consensus_track_file)
        except Exception as e:
            logger.error("Error reading %s: %s", consensus_track_file, str(e))
            logger.info("Backing up and recreating consensus track file")
            # Create backup before removing
            backup_file = consensus_track_file.with_suffix('.bed.bak')
            consensus_track_file.rename(backup_file)
            recompute_consensus_track = True
    else:
        if get_p_vals:
            recompute_consensus_track = True
    
    # If both files exist and are valid, and recompute is False, proceed to thresholding
    if not recompute_prior_localizations and not recompute_consensus_track:
        logger.info("Existing files are valid. Skipping recomputation.")
        
        if get_p_vals and alpha is not None:
            # Use the existing consensus_track.bed to generate significant regions
            sorted_df = track_df.copy()
            sorted_df['is_significant'] = sorted_df['p_value'] <= alpha
            sorted_df['block_change'] = sorted_df['is_significant'] != sorted_df['is_significant'].shift(1, fill_value=False)
            sorted_df['block_id'] = sorted_df['block_change'].cumsum()
            
            significant_blocks = sorted_df[sorted_df['is_significant']].groupby('block_id')
            
            if len(significant_blocks) > 0:
                bed_regions = significant_blocks.agg({
                    'chromosome': 'first',
                    'start': 'first',
                    'end': 'last'
                }).reset_index(drop=True)
                
                logger.info("Found %d significant consensus regions", len(bed_regions))
                
                bed_regions.to_csv(consensus_regions_file, sep='\t', header=False, index=False)
                logger.info("Consensus regions saved to %s", consensus_regions_file)
            else:
                logger.info("No significant consensus regions found")
        
        return  # Exit if no recomputation needed
    
    # If recomputation is needed, proceed
    logger.info("Recomputing necessary outputs...")
    
    # Process column names
    if table_columns:
        names = table_columns.split(',')
    location_column = 'location' if table_columns and 'location' in names else 0
    
    try:
        # Load data
        selected_loci_df, basis_identifiers_df, num_basis_identifiers = load_basis_identifiers(
            location, sep, table_columns, location_column
        )
        
        # Create prior localization regions
        prior_localization_df, atoms_df = create_prior_localizations(selected_loci_df, num_basis_identifiers)
        
        # Save prior localization regions
        prior_localization_df.to_csv(prior_localization_regions_file, sep='\t', header=False, index=False)
        logger.info("Prior localization regions saved to %s", prior_localization_regions_file)
        
        # Exit if p-values are not requested
        if not get_p_vals:
            logger.info("P-values not requested. Exiting after generating prior localization regions.")
            return
        
        # Compute selectivity
        atoms_df, identifier_coverage_sums, prior_localizations, selectivity_values, n_unique_instances = compute_identifier_selectivity(
            atoms_df, num_basis_identifiers
        )
        
        # Update basis identifiers with selectivity
        basis_identifiers_df['selectivity'] = selectivity_values[:-1]  # Exclude union
        logger.info("Added selectivity values to basis identifiers")
        
        # Ensure no zero selectivities
        selectivities = np.clip(selectivity_values, 1e-10, None)
        
        # Extract identifier vectors as a matrix for vectorized operations
        identifier_matrix_with_union = np.vstack(atoms_df['identifier_vector_with_union'].values)
        
        # Calculate consensus statistics using vectorized operations
        logger.info("Computing consensus statistics")
        consensus_scores = []
        for row in identifier_matrix_with_union:
            sel_score, omit_score = compute_consensus_statistics(row, selectivities)
            consensus_scores.append((sel_score, omit_score))
        
        # Add consensus scores to atoms_df
        consensus_scores_array = np.array(consensus_scores)
        atoms_df['consensus_selection_score'] = consensus_scores_array[:, 0]
        atoms_df['consensus_omission_score'] = consensus_scores_array[:, 1]
        
        # Simulate null distributions
        sorted_selection_stats, selection_ecdf, sorted_omission_stats, omission_ecdf = simulate_independence_null_distribution(
            selectivities, n_simulations
        )
        
        # Compute p-values
        logger.info("Computing p-values")
        observed_consensus_selection_statistics = atoms_df['consensus_selection_score'].values
        observed_consensus_omission_statistics = atoms_df['consensus_omission_score'].values
        
        atoms_df['p_selection'] = [
            compute_p_value(x, sorted_selection_stats, selection_ecdf) 
            for x in observed_consensus_selection_statistics
        ]
        
        atoms_df['p_omission'] = [
            compute_p_value(x, sorted_omission_stats, omission_ecdf) 
            for x in observed_consensus_omission_statistics
        ]
        
        # Build the consensus track
        logger.info("Building consensus track")
        track_df = atoms_df.loc[:, [
            'chromosome', 'position', 'next_position',
            'consensus_selection_score', 'p_selection',
            'consensus_omission_score', 'p_omission'
        ]].copy()
        
        track_df = track_df.rename(columns={'position': 'start', 'next_position': 'end'})
        track_df = track_df[track_df['end'] > track_df['start']]  # ensure valid intervals
        
        # Save consensus track
        track_df.to_csv(consensus_track_file, sep='\t', header=True, index=False)
        logger.info("Consensus track with p-values saved to %s", consensus_track_file)
        
        # Apply multiple testing correction if alpha is provided
        if alpha is not None:
            atoms_df, consensus_regions_sel, consensus_regions_omit = apply_multiple_testing_correction(
                atoms_df, n_unique_instances, alpha, mtc
            )
            
            # Save significant regions
            save_consensus_results(
                output_dir=output_path,
                prior_localization_regions_df=prior_localization_df,
                consensus_track_df=track_df,
                consensus_regions_sel=consensus_regions_sel,
                consensus_regions_omit=consensus_regions_omit,
                alpha=alpha,
                mtc=mtc
            )
        else:
            logger.info("Alpha threshold not provided. Only prior localization regions and consensus track were generated.")
            
    except Exception as e:
        logger.exception("Error during consensus computation: %s", str(e))
        raise


if __name__ == "__main__":
    construct_consensus_beds()