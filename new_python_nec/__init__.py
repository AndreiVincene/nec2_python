"""
Python Implementation of Numerical Electromagnetics Code (NEC-2)

This is a complete Python conversion of the NEC-2 Fortran code developed at
Lawrence Livermore Lab by G. Burke and A. Poggio, with modifications by Arie Voors.

The implementation includes:
- Method of Moments (MoM) electromagnetic analysis
- Wire and surface patch modeling
- Ground plane effects (perfect, finite, Sommerfeld)
- Network analysis and loading
- Far-field and near-field calculations
- Frequency sweeps and optimization support

Author: Python conversion of original NEC-2 Fortran code
Version: 1.0.0
"""

from .constants import *
from .data_structures import *
from .geometry import *
from .matrix_operations import *
from .field_calculations import *
from .ground_effects import *
from .network_analysis import *
from .input_output import *
from .main_engine import *

__version__ = "1.0.0"
__author__ = "Python NEC-2 Implementation"
__all__ = [
    "NECEngine",
    "GeometryData",
    "MatrixData", 
    "FieldData",
    "GroundData",
    "NetworkData",
    "InputData",
    "OutputData"
] 