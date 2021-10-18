# Enable importing local packages when imported into another module,
# while also using them as scripts

import sys
import os

cur_dir = os.path.join(os.path.dirname(__file__))
sys.path.append(cur_dir)
