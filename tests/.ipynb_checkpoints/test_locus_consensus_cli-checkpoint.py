#!/usr/bin/env python
"""
Tests for the consensus framework CLI functionality.

This module contains tests for the command-line interface of the consensus framework.
"""

import os
import tempfile
import shutil
from pathlib import Path
import pytest
from click.testing import CliRunner

# Import the CLI function to test
from consensus_framework.locus_consensus import construct_consensus_beds


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Clean up after the test
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_input_file():
    """Create a temporary sample input file."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.tsv') as f:
        f.write("location\tsep\ttable_columns\n")
        
        # Create tiny BED files for testing
        bed_files = []
        for i in range(3):
            bed_path = f.name.replace('.tsv', f'_caller{i}.bed')
            with open(bed_path, 'w') as bed:
                bed.write(f"chr1\t100\t200\n")
                bed.write(f"chr1\t300\t400\n")
            bed_files.append(bed_path)
            
            # Add the BED file to the input table
            f.write(f"{bed_path}\t\\t\tchromosome,start,end\n")
    
    yield f.name
    
    # Clean up
    os.unlink(f.name)
    for bed_file in bed_files:
        if os.path.exists(bed_file):
            os.unlink(bed_file)


def test_cli_basic_run(sample_input_file, temp_output_dir):
    """Test basic CLI functionality."""
    runner = CliRunner()
    result = runner.invoke(construct_consensus_beds, [
        '--location', sample_input_file,
        '--output-dir', temp_output_dir,
        '--skip-p-vals'  # Skip p-values to speed up the test
    ])
    
    # Check that the command ran successfully
    assert result.exit_code == 0
    
    # Check that the output file was created
    assert Path(temp_output_dir, 'prior_localization_regions.bed').exists()


def test_cli_with_p_values(sample_input_file, temp_output_dir):
    """Test CLI with p-value computation."""
    runner = CliRunner()
    result = runner.invoke(construct_consensus_beds, [
        '--location', sample_input_file,
        '--output-dir', temp_output_dir,
        '--get-p-vals',
        '--n-simulations', '100'  # Use fewer simulations for faster testing
    ])
    
    # Check that the command ran successfully
    assert result.exit_code == 0
    
    # Check that output files were created
    assert Path(temp_output_dir, 'prior_localization_regions.bed').exists()
    assert Path(temp_output_dir, 'consensus_track.bed').exists()


def test_cli_with_significance_threshold(sample_input_file, temp_output_dir):
    """Test CLI with significance threshold."""
    runner = CliRunner()
    result = runner.invoke(construct_consensus_beds, [
        '--location', sample_input_file,
        '--output-dir', temp_output_dir,
        '--get-p-vals',
        '--n-simulations', '100',  # Use fewer simulations for faster testing
        '--alpha', '0.05',
        '--mtc', 'bh'
    ])
    
    # Check that the command ran successfully
    assert result.exit_code == 0
    
    # Check that all output files were created
    assert Path(temp_output_dir, 'prior_localization_regions.bed').exists()
    assert Path(temp_output_dir, 'consensus_track.bed').exists()
    # Note: consensus regions may not exist if no regions are significant


def test_cli_verbose_flag(sample_input_file, temp_output_dir):
    """Test CLI with verbose flag."""
    runner = CliRunner()
    result = runner.invoke(construct_consensus_beds, [
        '--location', sample_input_file,
        '--output-dir', temp_output_dir,
        '--skip-p-vals',
        '--verbose'
    ])
    
    # Check that the command ran successfully
    assert result.exit_code == 0
    
    # Verbose output should contain DEBUG messages
    assert "DEBUG" in result.output


def test_cli_invalid_input():
    """Test CLI with invalid input file."""
    runner = CliRunner()
    result = runner.invoke(construct_consensus_beds, [
        '--location', 'nonexistent_file.tsv',
        '--output-dir', 'some_dir'
    ])
    
    # Command should fail with non-zero exit code
    assert result.exit_code != 0


if __name__ == "__main__":
    pytest.main()