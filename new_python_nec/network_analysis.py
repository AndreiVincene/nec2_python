"""
Network Analysis Module for Python NEC-2 Implementation

This module contains all network-related functions converted from the Fortran code,
including NETWK, COUPLE, and other network analysis subroutines.
"""

import numpy as np
from typing import Tuple, List, Optional, Complex
from .constants import *
from .data_structures import *

# ============================================================================
# NETWORK ANALYSIS (NETWK)
# ============================================================================

def netwk(cm: np.ndarray, cb: np.ndarray, cc: np.ndarray, cd: np.ndarray,
          ip: np.ndarray, cur: np.ndarray, nec_data: NECData) -> None:
    """
    Network analysis routine (NETWK from Fortran)
    
    Args:
        cm: Main matrix
        cb: B matrix for NGF
        cc: C matrix for NGF
        cd: D matrix for NGF
        ip: Pivot array
        cur: Current vector
        nec_data: Main NEC data structure
    """
    if nec_data.network.nonet == 0:
        # No network - solve directly
        solves(cm, ip, cur, nec_data.geometry.n + 2*nec_data.geometry.m, 1, 
               nec_data.geometry.np, nec_data.geometry.n, nec_data.geometry.mp, 
               nec_data.geometry.m, 13, 13)
        return
        
    # Network analysis with transmission lines and loads
    write_network_data(nec_data)
    
    # Initialize network parameters
    zpnorm = 0.0
    if nec_data.excitation.iped == 1:
        zpnorm = nec_data.excitation.zpnorm
        
    # Process each network
    for i in range(nec_data.network.nonet):
        ntyp = nec_data.network.ntyp[i]
        iseg1 = nec_data.network.iseg1[i]
        iseg2 = nec_data.network.iseg2[i]
        
        if ntyp == 1:
            # Transmission line
            process_transmission_line(i, iseg1, iseg2, cur, nec_data)
        elif ntyp == 2:
            # Admittance matrix
            process_admittance_matrix(i, iseg1, iseg2, cur, nec_data)
        elif ntyp == 3:
            # Impedance matrix
            process_impedance_matrix(i, iseg1, iseg2, cur, nec_data)
    
    # Solve the modified system
    if nec_data.network.icase > 0:
        # NGF case
        solve_ngf_network(cm, cb, cc, cd, ip, cur, nec_data)
    else:
        # Standard case
        solves(cm, ip, cur, nec_data.geometry.n + 2*nec_data.geometry.m, 1,
               nec_data.geometry.np, nec_data.geometry.n, nec_data.geometry.mp,
               nec_data.geometry.m, 13, 13)
    
    # Calculate input impedance if requested
    if nec_data.excitation.iped > 0:
        calculate_input_impedance(cur, nec_data)

def process_transmission_line(net_idx: int, iseg1: int, iseg2: int, 
                             cur: np.ndarray, nec_data: NECData) -> None:
    """
    Process transmission line network (NT card)
    
    Args:
        net_idx: Network index
        iseg1, iseg2: Connected segments
        cur: Current vector
        nec_data: Main NEC data structure
    """
    # Get transmission line parameters
    x11r = nec_data.network.x11r[net_idx]
    x11i = nec_data.network.x11i[net_idx]
    x12r = nec_data.network.x12r[net_idx]
    x12i = nec_data.network.x12i[net_idx]
    x22r = nec_data.network.x22r[net_idx]
    x22i = nec_data.network.x22i[net_idx]
    
    # Calculate transmission line length if not specified
    if x11i <= 0.0:
        x1 = nec_data.geometry.x[iseg1 - 1]
        y1 = nec_data.geometry.y[iseg1 - 1]
        z1 = nec_data.geometry.z[iseg1 - 1]
        x2 = nec_data.geometry.x[iseg2 - 1]
        y2 = nec_data.geometry.y[iseg2 - 1]
        z2 = nec_data.geometry.z[iseg2 - 1]
        
        length = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
        x11i = nec_data.geometry.wlam * length
    
    # Calculate transmission line parameters
    z0 = complex(x11r, x11i)
    length = x11i / nec_data.geometry.wlam
    beta = 2 * PI * length
    
    # Calculate ABCD parameters
    if abs(z0) > 0:
        cos_beta = np.cos(beta)
        sin_beta = np.sin(beta)
        
        a = cos_beta
        b = 1j * z0 * sin_beta
        c = 1j * sin_beta / z0
        d = cos_beta
        
        # Convert to Y parameters
        det = a * d - b * c
        if abs(det) > 1e-10:
            y11 = d / det
            y12 = -b / det
            y21 = -c / det
            y22 = a / det
        else:
            y11 = y12 = y21 = y22 = 0.0
    else:
        y11 = y12 = y21 = y22 = 0.0
    
    # Store Y parameters
    nec_data.network.x11r[net_idx] = np.real(y11)
    nec_data.network.x11i[net_idx] = np.imag(y11)
    nec_data.network.x12r[net_idx] = np.real(y12)
    nec_data.network.x12i[net_idx] = np.imag(y12)
    nec_data.network.x22r[net_idx] = np.real(y22)
    nec_data.network.x22i[net_idx] = np.imag(y22)

def process_admittance_matrix(net_idx: int, iseg1: int, iseg2: int,
                             cur: np.ndarray, nec_data: NECData) -> None:
    """
    Process admittance matrix network (TL card)
    
    Args:
        net_idx: Network index
        iseg1, iseg2: Connected segments
        cur: Current vector
        nec_data: Main NEC data structure
    """
    # Y parameters are already stored in the network data
    # No additional processing needed for admittance matrix
    pass

def process_impedance_matrix(net_idx: int, iseg1: int, iseg2: int,
                            cur: np.ndarray, nec_data: NECData) -> None:
    """
    Process impedance matrix network
    
    Args:
        net_idx: Network index
        iseg1, iseg2: Connected segments
        cur: Current vector
        nec_data: Main NEC data structure
    """
    # Convert Z parameters to Y parameters
    z11 = complex(nec_data.network.x11r[net_idx], nec_data.network.x11i[net_idx])
    z12 = complex(nec_data.network.x12r[net_idx], nec_data.network.x12i[net_idx])
    z22 = complex(nec_data.network.x22r[net_idx], nec_data.network.x22i[net_idx])
    
    det = z11 * z22 - z12 * z12
    if abs(det) > 1e-10:
        y11 = z22 / det
        y12 = -z12 / det
        y22 = z11 / det
    else:
        y11 = y12 = y22 = 0.0
    
    # Store Y parameters
    nec_data.network.x11r[net_idx] = np.real(y11)
    nec_data.network.x11i[net_idx] = np.imag(y11)
    nec_data.network.x12r[net_idx] = np.real(y12)
    nec_data.network.x12i[net_idx] = np.imag(y12)
    nec_data.network.x22r[net_idx] = np.real(y22)
    nec_data.network.x22i[net_idx] = np.imag(y22)

def solve_ngf_network(cm: np.ndarray, cb: np.ndarray, cc: np.ndarray, cd: np.ndarray,
                      ip: np.ndarray, cur: np.ndarray, nec_data: NECData) -> None:
    """
    Solve network for NGF case
    
    Args:
        cm: Main matrix
        cb: B matrix for NGF
        cc: C matrix for NGF
        cd: D matrix for NGF
        ip: Pivot array
        cur: Current vector
        nec_data: Main NEC data structure
    """
    # This is a simplified implementation
    # The full implementation would handle the complex NGF matrix operations
    
    # For now, use standard solver
    neq = nec_data.geometry.n + 2 * nec_data.geometry.m
    solves(cm, ip, cur, neq, 1, nec_data.geometry.np, nec_data.geometry.n,
           nec_data.geometry.mp, nec_data.geometry.m, 13, 13)

def calculate_input_impedance(cur: np.ndarray, nec_data: NECData) -> None:
    """
    Calculate input impedance
    
    Args:
        cur: Current vector
        nec_data: Main NEC data structure
    """
    if nec_data.excitation.nsant > 0:
        # Voltage source
        iseg = nec_data.excitation.isant[0]
        vs = nec_data.excitation.vsant[0]
        cur_seg = cur[iseg - 1]
        
        if abs(cur_seg) > 1e-20:
            zped = vs / cur_seg
        else:
            zped = complex(0.0, 0.0)
    elif nec_data.excitation.nvqd > 0:
        # Current source
        iseg = nec_data.excitation.ivqd[0]
        iq = nec_data.excitation.vqd[0]
        cur_seg = cur[iseg - 1]
        
        if abs(iq) > 1e-20:
            zped = cur_seg / iq
        else:
            zped = complex(0.0, 0.0)
    else:
        zped = complex(0.0, 0.0)
    
    nec_data.network.zped = zped

def write_network_data(nec_data: NECData) -> None:
    """Write network data to output"""
    print("\n" + "="*60)
    print("NETWORK DATA")
    print("="*60)
    
    for i in range(nec_data.network.nonet):
        ntyp = nec_data.network.ntyp[i]
        iseg1 = nec_data.network.iseg1[i]
        iseg2 = nec_data.network.iseg2[i]
        
        if ntyp == 1:
            print(f"Transmission Line {i+1}: Segments {iseg1} to {iseg2}")
        elif ntyp == 2:
            print(f"Admittance Matrix {i+1}: Segments {iseg1} to {iseg2}")
        elif ntyp == 3:
            print(f"Impedance Matrix {i+1}: Segments {iseg1} to {iseg2}")
            
        x11 = complex(nec_data.network.x11r[i], nec_data.network.x11i[i])
        x12 = complex(nec_data.network.x12r[i], nec_data.network.x12i[i])
        x22 = complex(nec_data.network.x22r[i], nec_data.network.x22i[i])
        
        print(f"  Y11 = {x11:.6f}")
        print(f"  Y12 = {x12:.6f}")
        print(f"  Y22 = {x22:.6f}")

# ============================================================================
# COUPLING ANALYSIS (COUPLE)
# ============================================================================

def couple(cur: np.ndarray, wlam: float, nec_data: NECData) -> None:
    """
    Compute maximum coupling between pairs of segments (COUPLE from Fortran)
    
    Args:
        cur: Current vector
        wlam: Wavelength
        nec_data: Main NEC data structure
    """
    if nec_data.network.nsant != 1 or nec_data.network.nvqd != 0:
        return
        
    # Get coupling parameters
    ncoupl = nec_data.network.ncoup
    if ncoupl == 0:
        return
        
    # Check if this is the right source segment
    j = isegno(nec_data.network.nctag[nec_data.network.icoup], 
               nec_data.network.ncseg[nec_data.network.icoup], nec_data)
    if j != nec_data.network.isant[0]:
        return
        
    # Increment coupling counter
    nec_data.network.icoup += 1
    
    # Calculate input impedance
    zin = nec_data.network.vsant[0]
    y11a = cur[j-1] * wlam / zin
    nec_data.network.y11a[nec_data.network.icoup - 1] = y11a
    
    # Calculate coupling coefficients
    l1 = (nec_data.network.icoup - 1) * (ncoupl - 1)
    for i in range(ncoupl):
        if i == nec_data.network.icoup - 1:
            continue
        k = isegno(nec_data.network.nctag[i], nec_data.network.ncseg[i], nec_data)
        l1 += 1
        nec_data.network.y12a[l1 - 1] = cur[k-1] * wlam / zin
    
    # If all couplings calculated, print results
    if nec_data.network.icoup >= ncoupl:
        print_coupling_results(nec_data)

def print_coupling_results(nec_data: NECData) -> None:
    """Print coupling analysis results"""
    print("\n" + "="*60)
    print("ISOLATION DATA")
    print("="*60)
    print("COUPLING BETWEEN SEGMENTS")
    print("TAG/SEG  NO.  TAG/SEG  NO.  COUPLING(DB)  LOAD IMPEDANCE(2ND SEG)  INPUT IMPEDANCE")
    print("                        REAL    IMAG.    REAL    IMAG.")
    print("-" * 80)
    
    ncoupl = nec_data.network.ncoup
    npm1 = ncoupl - 1
    
    for i in range(npm1):
        itt1 = nec_data.network.nctag[i]
        its1 = nec_data.network.ncseg[i]
        isg1 = isegno(itt1, its1, nec_data)
        
        for j in range(i + 1, ncoupl):
            itt2 = nec_data.network.nctag[j]
            its2 = nec_data.network.ncseg[j]
            isg2 = isegno(itt2, its2, nec_data)
            
            j1 = j + (i - 1) * npm1 - 1
            j2 = i + (j - 1) * npm1
            
            y11 = nec_data.network.y11a[i]
            y22 = nec_data.network.y11a[j]
            y12 = 0.5 * (nec_data.network.y12a[j1] + nec_data.network.y12a[j2])
            
            yin = y12 * y12
            dbc = abs(yin)
            c = dbc / (2.0 * np.real(y11) * np.real(y22) - np.real(yin))
            
            if c < 0.0 or c > 1.0:
                print(f"{itt1:4d} {its1:4d} {isg1:5d} {itt2:4d} {its2:4d} {isg2:5d} "
                      f"**ERROR** COUPLING IS NOT BETWEEN 0 AND 1. (= {c:.3e})")
                continue
                
            if c < 0.01:
                gmax = 0.5 * (c + 0.25 * c * c * c)
            else:
                gmax = (1.0 - np.sqrt(1.0 - c * c)) / c
                
            rho = gmax * np.conj(yin) / dbc
            yl = ((1.0 - rho) / (1.0 + rho) + 1.0) * np.real(y22) - y22
            zl = 1.0 / yl
            yin = y11 - yin / (y22 + yl)
            zin = 1.0 / yin
            dbc = db10(gmax)
            
            print(f"{itt1:4d} {its1:4d} {isg1:5d} {itt2:4d} {its2:4d} {isg2:5d} "
                  f"{dbc:9.3f} {np.real(zl):12.5e} {np.imag(zl):12.5e} "
                  f"{np.real(zin):12.5e} {np.imag(zin):12.5e}")

def isegno(itag: int, iseg: int, nec_data: NECData) -> int:
    """
    Find segment number from tag and segment
    
    Args:
        itag: Tag number
        iseg: Segment number
        nec_data: Main NEC data structure
        
    Returns:
        Segment index
    """
    if itag == 0:
        return iseg
    
    # Search for segment with matching tag
    for i in range(nec_data.geometry.n):
        if nec_data.geometry.itag[i] == itag:
            return i + 1
    
    return iseg

# ============================================================================
# MATRIX SOLVING (SOLVES)
# ============================================================================

def solves(a: np.ndarray, ip: np.ndarray, b: np.ndarray, nrow: int, nrhs: int,
           np: int, n: int, mp: int, m: int, iu1: int, iu2: int) -> None:
    """
    Solve matrix equation (SOLVES from Fortran)
    
    Args:
        a: Coefficient matrix
        ip: Pivot array
        b: Right-hand side vector
        nrow: Number of rows
        nrhs: Number of right-hand sides
        np: Number of segments in primary structure
        n: Total number of segments
        mp: Number of patches in primary structure
        m: Total number of patches
        iu1, iu2: Unit numbers for out-of-core storage
    """
    # This is a simplified implementation
    # The full implementation would handle out-of-core storage and symmetry
    
    neq = n + 2 * m
    npeq = np + 2 * mp
    
    # For in-core solution, use direct solver
    if nec_data.matrix.icase <= 2:
        # Direct solution
        for k in range(nrhs):
            ka = k * nrow
            solve_direct(a, ip, b[ka:ka+nrow], neq)
    else:
        # Out-of-core solution
        solve_out_of_core(a, ip, b, nrow, nrhs, neq, iu1, iu2)

def solve_direct(a: np.ndarray, ip: np.ndarray, b: np.ndarray, n: int) -> None:
    """
    Direct matrix solution using LU decomposition
    
    Args:
        a: Coefficient matrix (LU decomposed)
        ip: Pivot array
        b: Right-hand side vector
        n: Matrix size
    """
    # Forward substitution
    for i in range(n):
        ipi = ip[i] - 1
        if ipi != i:
            b[i], b[ipi] = b[ipi], b[i]
        
        if i < n - 1:
            for j in range(i + 1, n):
                b[j] -= a[j, i] * b[i]
    
    # Backward substitution
    for i in range(n - 1, -1, -1):
        if i < n - 1:
            for j in range(i + 1, n):
                b[i] -= a[i, j] * b[j]
        b[i] /= a[i, i]

def solve_out_of_core(a: np.ndarray, ip: np.ndarray, b: np.ndarray, 
                      nrow: int, nrhs: int, neq: int, iu1: int, iu2: int) -> None:
    """
    Out-of-core matrix solution
    
    Args:
        a: Coefficient matrix
        ip: Pivot array
        b: Right-hand side vector
        nrow: Number of rows
        nrhs: Number of right-hand sides
        neq: Number of equations
        iu1, iu2: Unit numbers for storage
    """
    # This is a placeholder for out-of-core solution
    # The full implementation would handle file I/O and block operations
    
    # For now, use direct solution
    for k in range(nrhs):
        ka = k * nrow
        solve_direct(a, ip, b[ka:ka+nrow], neq)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def db10(x: float) -> float:
    """
    Convert magnitude to dB (10*log10)
    
    Args:
        x: Magnitude value
        
    Returns:
        Value in dB
    """
    if x < 1e-20:
        return -999.99
    return 10.0 * np.log10(x)

def db20(x: float) -> float:
    """
    Convert magnitude to dB (20*log10)
    
    Args:
        x: Magnitude value
        
    Returns:
        Value in dB
    """
    if x < 1e-20:
        return -999.99
    return 20.0 * np.log10(x)

def calculate_power_budget(cur: np.ndarray, nec_data: NECData) -> None:
    """
    Calculate power budget
    
    Args:
        cur: Current vector
        nec_data: Main NEC data structure
    """
    # Calculate input power
    if nec_data.excitation.nsant > 0:
        iseg = nec_data.excitation.isant[0]
        vs = nec_data.excitation.vsant[0]
        cur_seg = cur[iseg - 1]
        pin = 0.5 * np.real(vs * np.conj(cur_seg))
    else:
        pin = 0.0
    
    # Calculate radiated power
    prad = 0.0
    for i in range(nec_data.geometry.n):
        cur_seg = cur[i]
        prad += 0.5 * abs(cur_seg)**2 * nec_data.geometry.si[i]
    
    # Calculate structure loss
    ploss = 0.0
    for i in range(nec_data.geometry.n):
        if abs(np.real(nec_data.network.zarray[i])) > 1e-20:
            cur_seg = cur[i]
            ploss += 0.5 * abs(cur_seg)**2 * np.real(nec_data.network.zarray[i]) * nec_data.geometry.si[i]
    
    # Calculate network loss
    pnls = 0.0  # Placeholder for network loss calculation
    
    # Calculate efficiency
    if pin > 0:
        efficiency = 100.0 * (pin - ploss - pnls) / pin
    else:
        efficiency = 0.0
    
    # Store results
    nec_data.network.pin = pin
    nec_data.network.prad = prad
    nec_data.network.ploss = ploss
    nec_data.network.pnls = pnls
    nec_data.network.efficiency = efficiency
    
    # Print results
    print("\n" + "="*60)
    print("POWER BUDGET")
    print("="*60)
    print(f"INPUT POWER     = {pin:.3e} WATTS")
    print(f"RADIATED POWER  = {prad:.3e} WATTS")
    print(f"STRUCTURE LOSS  = {ploss:.3e} WATTS")
    print(f"NETWORK LOSS    = {pnls:.3e} WATTS")
    print(f"EFFICIENCY      = {efficiency:.2f} PERCENT") 