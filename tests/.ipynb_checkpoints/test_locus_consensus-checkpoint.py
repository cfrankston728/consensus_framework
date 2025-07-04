#!/usr/bin/env python
"""
Unit tests for the consensus framework core functions.

This module contains unit tests for the key functions in the consensus framework.
Run with pytest: pytest test_core.py
"""

import numpy as np
import pandas as pd
import pytest
from pathlib import Path
import tempfile
import os

# Import functions to test
from core import (
    fast_indices_to_binary_vector,
    is_all_zeros_identifier_vector_vectorized,
    compute_consensus_statistics,
    compute_p_value,
    simulate_independence_null_distribution
)


class TestUtilityFunctions:
    """Tests for utility functions in the core module."""
    
    def test_fast_indices_to_binary_vector(self):
        """Test fast_indices_to_binary_vector creates correct binary vectors."""
        # Test with empty indices
        result = fast_indices_to_binary_vector([], 5)
        expected = np.zeros(5, dtype=int)
        np.testing.assert_array_equal(result, expected)
        
        # Test with single index
        result = fast_indices_to_binary_vector([2], 5)
        expected = np.array([0, 0, 1, 0, 0], dtype=int)
        np.testing.assert_array_equal(result, expected)
        
        # Test with multiple indices
        result = fast_indices_to_binary_vector([0, 3, 4], 5)
        expected = np.array([1, 0, 0, 1, 1], dtype=int)
        np.testing.assert_array_equal(result, expected)
    
    def test_is_all_zeros_identifier_vector_vectorized(self):
        """Test is_all_zeros_identifier_vector_vectorized identifies zero vectors correctly."""
        # Create test matrix
        matrix = np.array([
            [0, 0, 0],  # All zeros
            [1, 0, 0],  # Not all zeros
            [0, 0, 0],  # All zeros
            [0, 1, 0],  # Not all zeros
        ])
        
        result = is_all_zeros_identifier_vector_vectorized(matrix)
        expected = np.array([True, False, True, False])
        np.testing.assert_array_equal(result, expected)
        
        # Test with fractional values
        matrix = np.array([
            [0.1, 0.2, 0.3],  # All < 0.5
            [0.6, 0.2, 0.1],  # Not all < 0.5
            [0.4, 0.4, 0.4],  # All < 0.5
            [0.1, 0.5, 0.1],  # Not all < 0.5 (equal to 0.5)
        ])
        
        result = is_all_zeros_identifier_vector_vectorized(matrix)
        expected = np.array([True, False, True, False])
        np.testing.assert_array_equal(result, expected)


class TestStatisticalFunctions:
    """Tests for statistical functions in the core module."""
    
    def test_compute_consensus_statistics(self):
        """Test compute_consensus_statistics calculates correct scores."""
        # Setup test data
        identifier_vector = np.array([1, 0, 1, 0, 1, 0])  # 3 selections, 3 omissions
        selectivities = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])  # Last is union
        
        # Calculate expected results
        expected_selection_score = -(np.log(0.2) + np.log(0.4) + np.log(0.6))
        expected_omission_score = -(np.log(1-0.3) + np.log(1-0.5) + np.log(1-0.7))
        
        # Call function
        sel_score, omit_score = compute_consensus_statistics(identifier_vector, selectivities)
        
        # Check results
        assert np.isclose(sel_score, expected_selection_score)
        assert np.isclose(omit_score, expected_omission_score)
    
    def test_compute_p_value(self):
        """Test compute_p_value calculates correct p-values."""
        # Setup test data
        sorted_stats = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        ecdf = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
        
        # Test with observation in the range
        assert compute_p_value(2.5, sorted_stats, ecdf) == 0.6  # 1 - 0.4
        assert compute_p_value(4.5, sorted_stats, ecdf) == 0.2  # 1 - 0.8
        
        # Test with observation below range
        assert compute_p_value(0.5, sorted_stats, ecdf) == 1.0  # 1 - 0
        
        # Test with observation above range
        assert compute_p_value(5.5, sorted_stats, ecdf) == 0.0  # 1 - 1.0
    
    def test_simulate_independence_null_distribution(self):
        """Test simulate_independence_null_distribution generates valid distributions."""
        # Setup test data
        selectivity_values = np.array([0.2, 0.3, 0.4, 0.5])  # Last is union
        n_simulations = 100
        
        # Call function
        sorted_selection_stats, selection_ecdf, sorted_omission_stats, omission_ecdf = simulate_independence_null_distribution(
            selectivity_values, n_simulations
        )
        
        # Check dimensions
        assert len(sorted_selection_stats) == n_simulations
        assert len(selection_ecdf) == n_simulations
        assert len(sorted_omission_stats) == n_simulations
        assert len(omission_ecdf) == n_simulations
        
        # Check values are sorted
        assert np.all(np.diff(sorted_selection_stats) >= 0)
        assert np.all(np.diff(sorted_omission_stats) >= 0)
        
        # Check ECDF values are in [0, 1]
        assert np.all((selection_ecdf >= 0) & (selection_ecdf <= 1))
        assert np.all((omission_ecdf >= 0) & (omission_ecdf <= 1))
        
        # Check ECDF is monotonically increasing
        assert np.all(np.diff(selection_ecdf) >= 0)
        assert np.all(np.diff(omission_ecdf) >= 0)


if __name__ == "__main__":
    pytest.main()