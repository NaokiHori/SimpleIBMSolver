import os
import sys

sys.path.append(os.path.abspath("./conf_params"))

project = "Simple IBM Solver"
author = "Naoki Hori"
copyright = f"2019-2023, {author}"

extensions = [
        "sphinx.ext.mathjax",
]

from mathjax_params import mathjax_path
from mathjax_params import mathjax3_config

