"""Output TSVs of Wikipedia views and PubMed citations for genes
"""

import argparse

from citations.citations import Citations
from views.views import Views

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
        self.fetch_hints()

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
