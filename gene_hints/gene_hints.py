"""Output TSVs gene popularity by Wikipedia views and PubMed citations
"""

import argparse

from citations.citations import Citations
from views.views import Views
from merge_hints import merge_hints

class GeneHints():
    def __init__(self, num_days, excludes=None):
        self.num_days = num_days
        self.excludes = excludes

    def fetch_hints(self):
        excludes = self.excludes
        if not excludes or "views" not in excludes:
            Views().run()
        if not excludes or "citations" not in excludes:
            Citations().run(self.num_days)

    def run(self):
        """Output TSVs gene popularity by Wikipedia views and PubMed citations
        """
        self.fetch_hints()

        print("\n")

        stem = "data/homo-sapiens-"
        cite_path = f"{stem}pubmed-citations.tsv"
        view_path = f"{stem}wikipedia-views.tsv"
        hint_path = f"{stem}gene-hints.tsv"
        merge_hints(cite_path, view_path, hint_path)

# Command-line handler
if __name__ == "__main__":
    usage = """
    python3 gene_hints/gene_hints.py --num-days 365
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=usage
    )
    parser.add_argument(
        "--num-days",
        type=int,
        help="Number of days to analyze",
        default=180
    )
    parser.add_argument(
        "--excludes",
        nargs="*",
        help="Data types to exclude",
        choices=["views", "citations"]
    )
    args = parser.parse_args()
    num_days = args.num_days
    excludes = args.excludes

    GeneHints(num_days, excludes).run()
