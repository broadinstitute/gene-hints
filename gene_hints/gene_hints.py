"""Output TSVs gene popularity by Wikipedia views and PubMed citations
"""

import argparse

from citations.citations import Citations
from views.views import Views
from merge_hints import merge_hints

class GeneHints():
    def __init__(self, days, sort_by="count", only=None, debug=0):
        self.days = days
        self.sort_by = sort_by
        self.only = only
        self.debug = debug

        if self.days == 180 and self.debug > 0:
            self.days = 2

    def call_subpipelines(self):
        only = self.only
        cache = self.debug
        if not only or "views" in only:
            Views(cache).run(self.sort_by, self.debug)
        if not only or "citations" in only:
            Citations(cache).run(self.days, self.sort_by)

    def run(self):
        """Output TSVs gene popularity by Wikipedia views and PubMed citations
        """
        self.call_subpipelines()

        print("\n")

        stem = "data/homo-sapiens-"
        cite_path = f"{stem}pubmed-citations.tsv"
        view_path = f"{stem}wikipedia-views.tsv"
        hint_path = f"{stem}gene-hints.tsv"
        merge_hints(cite_path, view_path, hint_path)

# Command-line handler
if __name__ == "__main__":
    usage = """
    python3 gene_hints/gene_hints.py --days 365
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=usage
    )
    parser.add_argument(
        "--days",
        type=int,
        help="Number of days to analyze.  (default: %(default)i)",
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
        help=(
            "Metric by which to sort PubMed citations.  (default: %(default)s)"
        ),
        choices=["count", "delta", "rank", "rank_delta"],
        default="count"
    )
    parser.add_argument(
        "--debug",
        help=(
            "Get fast but incomplete data.  Useful to develop.  Levels:" +
                "0: use default `days`, don't cache.  " +
                "1: use `days 2`, cache download but not compute.  " +
                "2: like `debug 1`, and cache intermediate compute.  " +
                "(default: %(default)i)"
        ),
        type=int,
        choices=[0, 1, 2],
        default=0
    )
    args = parser.parse_args()
    days = args.days
    sort_by = args.sort_by
    only = args.only
    debug = args.debug

    GeneHints(days, sort_by, only, debug).run()
