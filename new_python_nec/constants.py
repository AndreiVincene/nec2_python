"""
Constants and Parameters for Python NEC-2 Implementation

This module contains all the constants, parameters, and configuration values
from the original Fortran NEC-2 code, converted to Python.
"""

import numpy as np
from typing import Dict, Any

# ============================================================================
# SYSTEM PARAMETERS (from nec2dpar.inc)
# ============================================================================

# Maximum segment counts for different executable versions
MAXSEG_500 = 500
MAXSEG_1K5 = 1500  
MAXSEG_3K0 = 3000
MAXSEG_5K0 = 5000
MAXSEG_8K0 = 8000
MAXSEG_11K = 11000

# Default configuration (500 segments)
MAXSEG = MAXSEG_500
MAXMAT = MAXSEG  # Max number of 'in-core' allocation
LOADMX = MAXSEG // 10  # Max number of LD cards
NSMAX = 64  # Max number of EX cards  
NETMX = 64  # Max number of segs connected to NT/TL
JMAX = 60   # Max segments connected to a single segment or junction

# Matrix storage
IRESRV = MAXMAT ** 2

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

# Speed of light (m/s)
CVEL = 299.8

# Mathematical constants
PI = np.pi
TA = PI / 180.0  # Degrees to radians conversion factor
TD = 180.0 / PI  # Radians to degrees conversion factor

# Complex unit
FJ = 1j

# ============================================================================
# CARD TYPE IDENTIFIERS (from ATST array)
# ============================================================================

CARD_TYPES = {
    'CE': 'COMMENT',
    'FR': 'FREQUENCY', 
    'LD': 'LOADING',
    'GN': 'GROUND',
    'EX': 'EXCITATION',
    'NT': 'NETWORK',
    'XQ': 'EXECUTE',
    'NE': 'NEAR_FIELD',
    'GD': 'GROUND_REPRESENTATION',
    'RP': 'RADIATION_PATTERN',
    'CM': 'COMMENT',
    'NX': 'NEXT_STRUCTURE',
    'EN': 'END',
    'TL': 'TRANSMISSION_LINE',
    'PT': 'PRINT_CONTROL',
    'KH': 'MATRIX_INTEGRATION_LIMIT',
    'NH': 'NEAR_FIELD_H',
    'PQ': 'CHARGE_PRINT_CONTROL',
    'EK': 'EXTENDED_THIN_WIRE_KERNEL',
    'WG': 'WRITE_GREEN_FUNCTION',
    'CP': 'COUPLING',
    'PL': 'PLOT_FLAGS'
}

# ============================================================================
# GEOMETRY CARD TYPES (from ATST array in DATAGN)
# ============================================================================

GEOMETRY_CARDS = {
    'GW': 'WIRE',
    'GX': 'WIRE_ARC', 
    'GR': 'WIRE_GRID',
    'GS': 'WIRE_SPHERE',
    'GE': 'WIRE_ELLIPSE',
    'GM': 'WIRE_MOVE',
    'SP': 'SURFACE_PATCH',
    'SM': 'SURFACE_MOVE',
    'GF': 'SURFACE_FILE',
    'GA': 'SURFACE_ARC',
    'SC': 'SURFACE_CONE',
    'GC': 'SURFACE_CYLINDER',
    'GH': 'HELIX'
}

# ============================================================================
# POLARIZATION TYPES
# ============================================================================

POLARIZATION_TYPES = {
    0: 'LINEAR',
    1: 'RIGHT', 
    2: 'LEFT'
}

# ============================================================================
# NETWORK TYPES
# ============================================================================

NETWORK_TYPES = {
    1: 'TRANSMISSION_LINE',
    2: 'ADMITTANCE_MATRIX',
    3: 'IMPEDANCE_MATRIX'
}

# ============================================================================
# GROUND TYPES
# ============================================================================

GROUND_TYPES = {
    0: 'FREE_SPACE',
    1: 'PERFECT_GROUND',
    2: 'FINITE_GROUND_SOMMERFELD'
}

# ============================================================================
# EXCITATION TYPES
# ============================================================================

EXCITATION_TYPES = {
    0: 'VOLTAGE_SOURCE',
    1: 'CURRENT_SOURCE',
    2: 'PLANE_WAVE',
    3: 'CURRENT_MOMENT',
    4: 'DIPOLE_MOMENT',
    5: 'MAGNETIC_FRILL'
}

# ============================================================================
# DEFAULT VALUES
# ============================================================================

# Default frequency and wavelength
DEFAULT_FREQ_MHZ = CVEL
DEFAULT_WAVELENGTH = CVEL / DEFAULT_FREQ_MHZ

# Default matrix integration limit
DEFAULT_RKH = 1.0

# Default normalization factor
NORMF = 200

# ============================================================================
# NUMERICAL PARAMETERS
# ============================================================================

# Convergence criteria
CONVERGENCE_TOLERANCE = 1e-4
MAX_ITERATIONS = 20

# Integration parameters
ROMBERG_TOLERANCE = 1e-4
MAX_ROMBERG_STEPS = 131072
ROMBERG_NTS = 4

# Sommerfeld integration
SOMMERFELD_CRITICAL = 1e-4
SOMMERFELD_MAX_H = 20

# ============================================================================
# FILE UNITS (from Fortran)
# ============================================================================

FILE_UNITS = {
    'INPUT': 2,
    'OUTPUT': 3, 
    'PLOT': 8,
    'MATRIX': 11,
    'MATRIX_C': 12,
    'MATRIX_B': 14,
    'MATRIX_D': 15,
    'SOM2D': 21
}

# ============================================================================
# ERROR MESSAGES
# ============================================================================

ERROR_MESSAGES = {
    'INVALID_CARD': 'INCORRECT LABEL FOR A COMMENT CARD',
    'FAULTY_CARD': 'FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION',
    'LOAD_EXCEED': 'NUMBER OF LOADING CARDS EXCEEDS STORAGE ALLOTTED',
    'EXCITATION_EXCEED': 'NUMBER OF EXCITATION CARDS EXCEEDS STORAGE ALLOTTED',
    'NETWORK_EXCEED': 'NUMBER OF NETWORK CARDS EXCEEDS STORAGE ALLOTTED',
    'NEAR_FIELD_MULTI': 'WHEN MULTIPLE FREQUENCIES ARE REQUESTED, ONLY ONE NEAR FIELD CARD CAN BE USED',
    'NGF_NOT_ALLOWED': 'CARD IS NOT ALLOWED WITH N.G.F.',
    'NGF_IN_USE': 'N.G.F. IN USE. CANNOT WRITE NEW N.G.F.',
    'COUPLING_EXCEED': 'NUMBER OF SEGMENTS IN COUPLING CALCULATION (CP) EXCEEDS LIMIT',
    'GROUND_RADIAL': 'RADIAL WIRE G. S. APPROXIMATION MAY NOT BE USED WITH SOMMERFELD GROUND OPTION',
    'GROUND_PARAM_ERROR': 'ERROR IN GROUND PARAMETERS',
    'SOM2D_FILE_ERROR': 'ERROR OPENING SOMMERFELD GROUND FILE - SOM2D.NEC'
}

# ============================================================================
# FORMAT STRINGS (from Fortran FORMAT statements)
# ============================================================================

FORMAT_STRINGS = {
    'HEADER': '*********************************************\nNUMERICAL ELECTROMAGNETICS CODE (NEC-2D)\n*********************************************',
    'FREQUENCY': 'FREQUENCY={:.4e} MHZ\nWAVELENGTH={:.4e} METERS',
    'MATRIX_TIMING': 'FILL={:.3f} SEC., FACTOR={:.3f} SEC.',
    'POWER_BUDGET': 'INPUT POWER={:.4e} WATTS\nRADIATED POWER={:.4e} WATTS\nSTRUCTURE LOSS={:.4e} WATTS\nNETWORK LOSS={:.4e} WATTS\nEFFICIENCY={:.2f} PERCENT',
    'CURRENT_HEADER': 'SEG. TAG COORD. OF SEG. CENTER SEG. - - - CURRENT (AMPS) - - -\nNO. NO. X Y Z LENGTH REAL IMAG. MAG. PHASE',
    'IMPEDANCE_HEADER': 'FREQ. - UNNORMALIZED IMPEDANCE - - NORMALIZED IMPEDANCE -\nMHZ RESISTANCE REACTANCE MAGNITUDE PHASE RESISTANCE REACTANCE MAGNITUDE PHASE'
}

# ============================================================================
# CONFIGURATION FUNCTIONS
# ============================================================================

def set_max_segments(max_seg: int) -> None:
    """Set the maximum number of segments for the simulation."""
    global MAXSEG, MAXMAT, LOADMX, IRESRV
    MAXSEG = max_seg
    MAXMAT = max_seg
    LOADMX = max_seg // 10
    IRESRV = MAXMAT ** 2

def get_configuration() -> Dict[str, Any]:
    """Get current configuration parameters."""
    return {
        'MAXSEG': MAXSEG,
        'MAXMAT': MAXMAT, 
        'LOADMX': LOADMX,
        'NSMAX': NSMAX,
        'NETMX': NETMX,
        'JMAX': JMAX,
        'IRESRV': IRESRV,
        'CVEL': CVEL,
        'PI': PI,
        'FJ': FJ
    }

def print_configuration() -> None:
    """Print current configuration parameters."""
    config = get_configuration()
    print("NEC-2 Python Configuration:")
    print("=" * 40)
    for key, value in config.items():
        print(f"{key:>10}: {value}")
    print("=" * 40) 