"""Tests for Wikipedia view counts and derived metrics

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

from gene_hints.views.views import Views

def test_save_to_file(tmpdir):

    # About tmpdir:
    # https://docs.pytest.org/en/6.2.x/tmpdir.html#the-tmpdir-fixture
    views = Views(output_dir=tmpdir)

    # Simulate view counts for three genes
    day_ago = {"INS": 900, "OXT": 30, "FOO": 10}
    two_days_ago = {"INS": 901, "OXT": 25, "FOO": 60}

    # Save counts (and derived metrics) for yesterday and day before
    views.save_to_file(day_ago, 0)
    views.save_to_file(two_days_ago, 1)

    # Parse output
    output_path = tmpdir + "homo-sapiens-wikipedia-views.tsv"
    with open(output_path) as f:
        lines = [line.split("\t") for line in f.readlines()]

    # Verify things output as expected
    assert len(lines) == 4
    assert len(lines[0]) == 5 # Expect five columns
    expected_lines = [
        ['# gene', 'views', 'view_delta', 'view_rank', 'view_rank_delta\n'],
        ['INS', '901', '1', '1', '0\n'],
        ['FOO', '60', '50', '2', '-1\n'],
        ['OXT', '25', '-5', '3', '1\n']
    ]
    assert lines == expected_lines


