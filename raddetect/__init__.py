__version__ = "0.1.3"

# Explicitly import the submodules (using the 'as' alias to satisfy Ruff F401)
from . import blumchen as blumchen
from . import cryoradon as cryoradon
from . import models as models
from . import monalpha as monalpha

# Define __all__ to control what 'from raddetect import *' does
__all__ = [
    "blumchen",
    "cryoradon",
    "monalpha",
    "SinglePeak",
    "DoublePeak",
    "TriplePeak",
    "AccumulationModel",
    "ExpansionModel",
    "linear_model",
    "res_model",
    "calculate_fractions",
]
