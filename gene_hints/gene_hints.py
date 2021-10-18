"""Output TSVs gene popularity by Wikipedia views and PubMed citations
"""

import argparse

from citations.citations import Citations
from views.views import Views
from merge_hints import merge_hints

class GeneHints():
    def __init__(self, num_days, sort_by="count", only=None):
        self.num_days = num_days
        self.sort_by = sort_by
        self.only = only

    def fetch_hints(self):
        only = self.only
        if not only or "views" not in only:
            Views().run()
        if not only or "citations" not in only:
            Citations().run(self.num_days, self.sort_by)

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
        "--only",
        nargs="*",
        help="Data types to include",
        choices=["views", "citations"]
    )
    parser.add_argument(
        "--sort-by",
        help="Metric by which to sort PubMed citations",
        choices=["count", "delta", "rank", "rank_delta"],
        default="count"
    )
    args = parser.parse_args()
    num_days = args.num_days
    sort_by = args.sort_by
    only = args.only

    GeneHints(num_days, sort_by, only).run()
