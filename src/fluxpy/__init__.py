__author__ = "Haris Zafeiropoulos"
__version__ = "0.0.1"

import os
from pathlib import Path

from fluxpy import utils
from fluxpy import stats
from fluxpy import illustrations

from .auto_decorator import auto_decorate

auto_decorate()

__all__ = [
    'illustrations',
    'mappings',
    'stats',
    'utils',
]

root  = Path(os.path.realpath(__file__)).parent
check = root / "data/MSEED_COMPOUNDS.json"

if not os.path.isfile(check):
    print("Downloading ModelSEED compound and reaction files")
    os.system(f"python {root}/data/get_modelseed_biochem.py") 
 


