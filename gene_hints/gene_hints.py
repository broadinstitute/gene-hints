"""Output TSVs of Wikipedia views and PubMed citations for genes
"""

import argparse

from citations.citations import Citations
from views.views import Views

class GeneHints():
    def __init__(
        self,
        num_days=180
    ):
        self.num_days = num_days

    def run(self):
        views = Views().run()
        cites = Citations().run(self.num_days)

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
        help="Days in timeframe to analyze.  Defaults to 180."
    )
    args = parser.parse_args()
    num_days = args.num_days

    GeneHints().run()
