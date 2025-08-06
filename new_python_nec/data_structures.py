"""
Data Structures for Python NEC-2 Implementation

This module defines all the data structures and classes that correspond to
the COMMON blocks in the original Fortran NEC-2 code.
"""

import numpy as np
from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
from .constants import *

# ============================================================================
# MAIN DATA STRUCTURE (corresponds to COMMON /DATA/)
# ============================================================================

@dataclass
class GeometryData:
    """Main geometry data structure corresponding to COMMON /DATA/"""
    
    # Segment coordinates and dimensions
    x: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    y: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    z: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    si: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))  # Segment lengths
    bi: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))  # Segment radii
    
    # Segment direction cosines
    alp: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))  # Alpha
    bet: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))  # Beta
    
    # Wavelength
    wlam: float = DEFAULT_WAVELENGTH
    
    # Connection data
    icon1: np.ndarray = field(default_factory=lambda: np.zeros(2*MAXSEG, dtype=np.int32))
    icon2: np.ndarray = field(default_factory=lambda: np.zeros(2*MAXSEG, dtype=np.int32))
    itag: np.ndarray = field(default_factory=lambda: np.zeros(2*MAXSEG, dtype=np.int32))
    iconx: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.int32))
    
    # Segment counts
    ld: int = MAXSEG
    n1: int = 0
    n2: int = 1
    n: int = 0   # Total wire segments
    np: int = 0  # Wire segments in symmetric cell
    m1: int = 0
    m2: int = 1
    m: int = 0   # Total surface patches
    mp: int = 0  # Surface patches in symmetric cell
    
    # Symmetry flag
    ipsym: int = 0
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.x) != MAXSEG:
            self.x = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.y) != MAXSEG:
            self.y = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.z) != MAXSEG:
            self.z = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.si) != MAXSEG:
            self.si = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.bi) != MAXSEG:
            self.bi = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.alp) != MAXSEG:
            self.alp = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.bet) != MAXSEG:
            self.bet = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.icon1) != 2*MAXSEG:
            self.icon1 = np.zeros(2*MAXSEG, dtype=np.int32)
        if len(self.icon2) != 2*MAXSEG:
            self.icon2 = np.zeros(2*MAXSEG, dtype=np.int32)
        if len(self.itag) != 2*MAXSEG:
            self.itag = np.zeros(2*MAXSEG, dtype=np.int32)
        if len(self.iconx) != MAXSEG:
            self.iconx = np.zeros(MAXSEG, dtype=np.int32)

# ============================================================================
# MATRIX DATA STRUCTURE (corresponds to COMMON /CMB/ and /MATPAR/)
# ============================================================================

@dataclass
class MatrixData:
    """Matrix data structure corresponding to COMMON /CMB/ and /MATPAR/"""
    
    # Main interaction matrix
    cm: np.ndarray = field(default_factory=lambda: np.zeros(IRESRV, dtype=np.complex128))
    
    # Matrix parameters
    icase: int = 0
    nbloks: int = 0
    npblk: int = 0
    nlast: int = 0
    nblsym: int = 0
    npsym: int = 0
    nlsym: int = 0
    imat: int = 0
    icasx: int = 0
    nbbx: int = 0
    npbx: int = 0
    nlbx: int = 0
    nbbl: int = 0
    npbl: int = 0
    nlbl: int = 0
    
    def __post_init__(self):
        """Initialize matrix with proper size"""
        if len(self.cm) != IRESRV:
            self.cm = np.zeros(IRESRV, dtype=np.complex128)

# ============================================================================
# CURRENT DATA STRUCTURE (corresponds to COMMON /CRNT/)
# ============================================================================

@dataclass
class CurrentData:
    """Current data structure corresponding to COMMON /CRNT/"""
    
    # Current coefficients
    air: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    aii: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    bir: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    bii: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    cir: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    cii: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.float64))
    
    # Current vector
    cur: np.ndarray = field(default_factory=lambda: np.zeros(3*MAXSEG, dtype=np.complex128))
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.air) != MAXSEG:
            self.air = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.aii) != MAXSEG:
            self.aii = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.bir) != MAXSEG:
            self.bir = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.bii) != MAXSEG:
            self.bii = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.cir) != MAXSEG:
            self.cir = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.cii) != MAXSEG:
            self.cii = np.zeros(MAXSEG, dtype=np.float64)
        if len(self.cur) != 3*MAXSEG:
            self.cur = np.zeros(3*MAXSEG, dtype=np.complex128)

# ============================================================================
# GROUND DATA STRUCTURE (corresponds to COMMON /GND/)
# ============================================================================

@dataclass
class GroundData:
    """Ground data structure corresponding to COMMON /GND/"""
    
    # Ground parameters
    zrati: complex = 1.0 + 0j
    zrati2: complex = 1.0 + 0j
    frati: complex = 1.0 + 0j
    t1: complex = 0.0 + 0j
    t2: complex = 0.0 + 0j
    cl: float = 0.0
    ch: float = 0.0
    scrwl: float = 0.0
    scrwr: float = 0.0
    nradl: int = 0
    ksymp: int = 1
    ifar: int = -1
    iperf: int = 0

# ============================================================================
# LOADING DATA STRUCTURE (corresponds to COMMON /ZLOAD/)
# ============================================================================

@dataclass
class LoadingData:
    """Loading data structure corresponding to COMMON /ZLOAD/"""
    
    # Load impedances
    zarray: np.ndarray = field(default_factory=lambda: np.zeros(MAXSEG, dtype=np.complex128))
    nload: int = 0
    nlodf: int = 0
    
    # Load parameters
    ldtyp: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.int32))
    ldtag: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.int32))
    ldtagf: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.int32))
    ldtagt: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.int32))
    zlr: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.float64))
    zli: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.float64))
    zlc: np.ndarray = field(default_factory=lambda: np.zeros(LOADMX, dtype=np.float64))
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.zarray) != MAXSEG:
            self.zarray = np.zeros(MAXSEG, dtype=np.complex128)
        if len(self.ldtyp) != LOADMX:
            self.ldtyp = np.zeros(LOADMX, dtype=np.int32)
        if len(self.ldtag) != LOADMX:
            self.ldtag = np.zeros(LOADMX, dtype=np.int32)
        if len(self.ldtagf) != LOADMX:
            self.ldtagf = np.zeros(LOADMX, dtype=np.int32)
        if len(self.ldtagt) != LOADMX:
            self.ldtagt = np.zeros(LOADMX, dtype=np.int32)
        if len(self.zlr) != LOADMX:
            self.zlr = np.zeros(LOADMX, dtype=np.float64)
        if len(self.zli) != LOADMX:
            self.zli = np.zeros(LOADMX, dtype=np.float64)
        if len(self.zlc) != LOADMX:
            self.zlc = np.zeros(LOADMX, dtype=np.float64)

# ============================================================================
# EXCITATION DATA STRUCTURE (corresponds to COMMON /VSORC/)
# ============================================================================

@dataclass
class ExcitationData:
    """Excitation data structure corresponding to COMMON /VSORC/"""
    
    # Voltage sources
    vqd: np.ndarray = field(default_factory=lambda: np.zeros(NSMAX, dtype=np.complex128))
    vsant: np.ndarray = field(default_factory=lambda: np.zeros(NSMAX, dtype=np.complex128))
    vqds: np.ndarray = field(default_factory=lambda: np.zeros(NSMAX, dtype=np.complex128))
    ivqd: np.ndarray = field(default_factory=lambda: np.zeros(NSMAX, dtype=np.int32))
    isant: np.ndarray = field(default_factory=lambda: np.zeros(NSMAX, dtype=np.int32))
    iqds: np.ndarray = field(default_factory=lambda: np.zeros(NSMAX, dtype=np.int32))
    nvqd: int = 0
    nsant: int = 0
    nqds: int = 0
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.vqd) != NSMAX:
            self.vqd = np.zeros(NSMAX, dtype=np.complex128)
        if len(self.vsant) != NSMAX:
            self.vsant = np.zeros(NSMAX, dtype=np.complex128)
        if len(self.vqds) != NSMAX:
            self.vqds = np.zeros(NSMAX, dtype=np.complex128)
        if len(self.ivqd) != NSMAX:
            self.ivqd = np.zeros(NSMAX, dtype=np.int32)
        if len(self.isant) != NSMAX:
            self.isant = np.zeros(NSMAX, dtype=np.int32)
        if len(self.iqds) != NSMAX:
            self.iqds = np.zeros(NSMAX, dtype=np.int32)

# ============================================================================
# NETWORK DATA STRUCTURE (corresponds to COMMON /NETCX/)
# ============================================================================

@dataclass
class NetworkData:
    """Network data structure corresponding to COMMON /NETCX/"""
    
    # Network parameters
    zped: complex = 0.0 + 0j
    pin: float = 0.0
    pnls: float = 0.0
    
    # Admittance matrix elements
    x11r: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.float64))
    x11i: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.float64))
    x12r: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.float64))
    x12i: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.float64))
    x22r: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.float64))
    x22i: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.float64))
    
    # Network types and connections
    ntyp: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.int32))
    iseg1: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.int32))
    iseg2: np.ndarray = field(default_factory=lambda: np.zeros(NETMX, dtype=np.int32))
    
    # Network counts
    neq: int = 0
    npeq: int = 0
    neq2: int = 0
    nonet: int = 0
    ntsol: int = 0
    nprint: int = 0
    masym: int = 0
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.x11r) != NETMX:
            self.x11r = np.zeros(NETMX, dtype=np.float64)
        if len(self.x11i) != NETMX:
            self.x11i = np.zeros(NETMX, dtype=np.float64)
        if len(self.x12r) != NETMX:
            self.x12r = np.zeros(NETMX, dtype=np.float64)
        if len(self.x12i) != NETMX:
            self.x12i = np.zeros(NETMX, dtype=np.float64)
        if len(self.x22r) != NETMX:
            self.x22r = np.zeros(NETMX, dtype=np.float64)
        if len(self.x22i) != NETMX:
            self.x22i = np.zeros(NETMX, dtype=np.float64)
        if len(self.ntyp) != NETMX:
            self.ntyp = np.zeros(NETMX, dtype=np.int32)
        if len(self.iseg1) != NETMX:
            self.iseg1 = np.zeros(NETMX, dtype=np.int32)
        if len(self.iseg2) != NETMX:
            self.iseg2 = np.zeros(NETMX, dtype=np.int32)

# ============================================================================
# SEGMENT JUNCTION DATA STRUCTURE (corresponds to COMMON /SEGJ/)
# ============================================================================

@dataclass
class SegmentJunctionData:
    """Segment junction data structure corresponding to COMMON /SEGJ/"""
    
    # Junction coefficients
    ax: np.ndarray = field(default_factory=lambda: np.zeros(JMAX, dtype=np.float64))
    bx: np.ndarray = field(default_factory=lambda: np.zeros(JMAX, dtype=np.float64))
    cx: np.ndarray = field(default_factory=lambda: np.zeros(JMAX, dtype=np.float64))
    jco: np.ndarray = field(default_factory=lambda: np.zeros(JMAX, dtype=np.int32))
    
    # Junction parameters
    jsno: int = 0
    iscon: np.ndarray = field(default_factory=lambda: np.zeros(50, dtype=np.int32))
    nscon: int = 0
    ipcon: np.ndarray = field(default_factory=lambda: np.zeros(10, dtype=np.int32))
    npcon: int = 0
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.ax) != JMAX:
            self.ax = np.zeros(JMAX, dtype=np.float64)
        if len(self.bx) != JMAX:
            self.bx = np.zeros(JMAX, dtype=np.float64)
        if len(self.cx) != JMAX:
            self.cx = np.zeros(JMAX, dtype=np.float64)
        if len(self.jco) != JMAX:
            self.jco = np.zeros(JMAX, dtype=np.int32)
        if len(self.iscon) != 50:
            self.iscon = np.zeros(50, dtype=np.int32)
        if len(self.ipcon) != 10:
            self.ipcon = np.zeros(10, dtype=np.int32)

# ============================================================================
# FIELD PATTERN DATA STRUCTURE (corresponds to COMMON /FPAT/)
# ============================================================================

@dataclass
class FieldPatternData:
    """Field pattern data structure corresponding to COMMON /FPAT/"""
    
    # Pattern angles
    thets: float = 0.0
    phis: float = 0.0
    dth: float = 1.0
    dph: float = 0.0
    
    # Field parameters
    rfld: float = 0.0
    gnor: float = 1.0
    clt: float = 0.0
    cht: float = 0.0
    epsr2: float = 1.0
    sig2: float = 0.0
    xpr6: float = 0.0
    pinr: float = 0.0
    pnlr: float = 0.0
    ploss: float = 0.0
    
    # Near field parameters
    xnr: float = 0.0
    ynr: float = 0.0
    znr: float = 0.0
    dxnr: float = 0.0
    dynr: float = 0.0
    dznr: float = 0.0
    nth: int = 91
    nph: int = 1
    ipd: int = 0
    iavp: int = 0
    inor: int = 0
    iax: int = 0
    ixtyp: int = 0
    near: int = -1
    nfeh: int = 0
    nrx: int = 0
    nry: int = 0
    nrz: int = 0

# ============================================================================
# COUPLING DATA STRUCTURE (corresponds to COMMON /YPARM/)
# ============================================================================

@dataclass
class CouplingData:
    """Coupling data structure corresponding to COMMON /YPARM/"""
    
    # Coupling parameters
    y11a: np.ndarray = field(default_factory=lambda: np.zeros(5, dtype=np.complex128))
    y12a: np.ndarray = field(default_factory=lambda: np.zeros(20, dtype=np.complex128))
    ncoup: int = 0
    icoup: int = 0
    nctag: np.ndarray = field(default_factory=lambda: np.zeros(5, dtype=np.int32))
    ncseg: np.ndarray = field(default_factory=lambda: np.zeros(5, dtype=np.int32))
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.y11a) != 5:
            self.y11a = np.zeros(5, dtype=np.complex128)
        if len(self.y12a) != 20:
            self.y12a = np.zeros(20, dtype=np.complex128)
        if len(self.nctag) != 5:
            self.nctag = np.zeros(5, dtype=np.int32)
        if len(self.ncseg) != 5:
            self.ncseg = np.zeros(5, dtype=np.int32)

# ============================================================================
# MAIN NEC DATA STRUCTURE
# ============================================================================

@dataclass
class NECData:
    """Main data structure containing all NEC-2 data"""
    
    # Core data structures
    geometry: GeometryData = field(default_factory=GeometryData)
    matrix: MatrixData = field(default_factory=MatrixData)
    current: CurrentData = field(default_factory=CurrentData)
    ground: GroundData = field(default_factory=GroundData)
    loading: LoadingData = field(default_factory=LoadingData)
    excitation: ExcitationData = field(default_factory=ExcitationData)
    network: NetworkData = field(default_factory=NetworkData)
    junction: SegmentJunctionData = field(default_factory=SegmentJunctionData)
    pattern: FieldPatternData = field(default_factory=FieldPatternData)
    coupling: CouplingData = field(default_factory=CouplingData)
    
    # Additional data (from other COMMON blocks)
    epsr: float = 1.0
    sig: float = 0.0
    scrwlt: float = 0.0
    scrwrt: float = 0.0
    fmhz: float = DEFAULT_FREQ_MHZ
    ip: np.ndarray = field(default_factory=lambda: np.zeros(2*MAXSEG, dtype=np.int32))
    kcom: int = 0
    
    # Comments
    com: np.ndarray = field(default_factory=lambda: np.zeros((19, 5), dtype='U4'))
    
    # Plot flags
    iplp1: int = 0
    iplp2: int = 0
    iplp3: int = 0
    iplp4: int = 0
    
    def __post_init__(self):
        """Initialize arrays with proper sizes"""
        if len(self.ip) != 2*MAXSEG:
            self.ip = np.zeros(2*MAXSEG, dtype=np.int32)
        if self.com.shape != (19, 5):
            self.com = np.zeros((19, 5), dtype='U4')
    
    def reset(self):
        """Reset all data to initial state"""
        self.geometry = GeometryData()
        self.matrix = MatrixData()
        self.current = CurrentData()
        self.ground = GroundData()
        self.loading = LoadingData()
        self.excitation = ExcitationData()
        self.network = NetworkData()
        self.junction = SegmentJunctionData()
        self.pattern = FieldPatternData()
        self.coupling = CouplingData()
        
        self.epsr = 1.0
        self.sig = 0.0
        self.scrwlt = 0.0
        self.scrwrt = 0.0
        self.fmhz = DEFAULT_FREQ_MHZ
        self.ip = np.zeros(2*MAXSEG, dtype=np.int32)
        self.kcom = 0
        self.com = np.zeros((19, 5), dtype='U4')
        self.iplp1 = 0
        self.iplp2 = 0
        self.iplp3 = 0
        self.iplp4 = 0 