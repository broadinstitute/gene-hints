"""Tests for PubMed citation counts and derived metrics

To run:
    $ pwd
    gene_hints
    $ cd tests
    $ pytest -s
"""

import sys

# Ensures `gene_hints` packages (and subpackages, like `views`) can be imported
# TODO: Find way to avoid this kludge
sys.path += ['..', '../gene_hints']

from gene_hints.citations.enrich_citations import rank_counts

def test_rank_counts():
    counts_by_key = {
        '59323': 7, '112400': 7, '293615': 7, '140914': 7, '361630': 5,
        '24318': 9, '310553': 7, '364594': 7, '24221': 7, '25504': 7,
        '24225': 8, '29705': 8, '29709': 6
    }

    ranks_by_key = rank_counts(counts_by_key)

    print('ranks_by_key', ranks_by_key)
    expected_ranks_by_key = {
        '59323': 4, '112400': 4, '293615': 4, '140914': 4, '361630': 13,
        '24318': 1, '310553': 4, '364594': 4, '24221': 4, '25504': 4,
        '24225': 2, '29705': 2, '29709': 12
    }

    assert ranks_by_key == expected_ranks_by_key

