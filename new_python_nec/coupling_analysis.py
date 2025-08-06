"""
Coupling Analysis Module for Python NEC-2 Implementation

This module contains the coupling analysis functions converted from the Fortran code,
including COUPLE and related subroutines for computing maximum coupling between segments.
"""

import numpy as np
from typing import List, Tuple, Optional, Complex
from .constants import *
from .data_structures import *

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
    # Check if coupling analysis is requested
    if nec_data.coupling.ncoup == 0:
        return
    
    # Check if we have the right source configuration
    if nec_data.excitation.nsant != 1 or nec_data.excitation.nvqd != 0:
        return
    
    # Get coupling parameters
    ncoupl = nec_data.coupling.ncoup
    icoup = nec_data.coupling.icoup
    
    # Check if this is the right source segment
    j = isegno(nec_data.coupling.nctag[icoup], nec_data.coupling.ncseg[icoup], nec_data)
    if j != nec_data.excitation.isant[0]:
        return
    
    # Increment coupling counter
    nec_data.coupling.icoup += 1
    
    # Calculate input impedance
    zin = nec_data.excitation.vsant[0]
    y11a = cur[j-1] * wlam / zin
    nec_data.coupling.y11a[nec_data.coupling.icoup - 1] = y11a
    
    # Calculate coupling coefficients
    l1 = (nec_data.coupling.icoup - 1) * (ncoupl - 1)
    for i in range(ncoupl):
        if i == nec_data.coupling.icoup - 1:
            continue
        k = isegno(nec_data.coupling.nctag[i], nec_data.coupling.ncseg[i], nec_data)
        l1 += 1
        if l1 < len(nec_data.coupling.y12a):
            nec_data.coupling.y12a[l1 - 1] = cur[k-1] * wlam / zin
    
    # If all couplings calculated, print results
    if nec_data.coupling.icoup >= ncoupl:
        print_coupling_results(nec_data)

def print_coupling_results(nec_data: NECData) -> None:
    """
    Print coupling analysis results
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\n" + "="*60)
    print("ISOLATION DATA")
    print("="*60)
    print("COUPLING BETWEEN SEGMENTS")
    print("TAG/SEG  NO.  TAG/SEG  NO.  COUPLING(DB)  LOAD IMPEDANCE(2ND SEG)  INPUT IMPEDANCE")
    print("                        REAL    IMAG.    REAL    IMAG.")
    print("-" * 80)
    
    ncoupl = nec_data.coupling.ncoup
    npm1 = ncoupl - 1
    
    for i in range(npm1):
        itt1 = nec_data.coupling.nctag[i]
        its1 = nec_data.coupling.ncseg[i]
        isg1 = isegno(itt1, its1, nec_data)
        
        for j in range(i + 1, ncoupl):
            itt2 = nec_data.coupling.nctag[j]
            its2 = nec_data.coupling.ncseg[j]
            isg2 = isegno(itt2, its2, nec_data)
            
            j1 = j + (i - 1) * npm1 - 1
            j2 = i + (j - 1) * npm1
            
            y11 = nec_data.coupling.y11a[i]
            y22 = nec_data.coupling.y11a[j]
            
            # Get coupling coefficient
            if j1 < len(nec_data.coupling.y12a) and j2 < len(nec_data.coupling.y12a):
                y12 = 0.5 * (nec_data.coupling.y12a[j1] + nec_data.coupling.y12a[j2])
            else:
                y12 = 0.0
            
            # Calculate coupling parameters
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
# COUPLING CALCULATION HELPERS
# ============================================================================

def calculate_coupling_coefficient(seg1: int, seg2: int, nec_data: NECData) -> complex:
    """
    Calculate coupling coefficient between two segments
    
    Args:
        seg1: First segment index
        seg2: Second segment index
        nec_data: Main NEC data structure
        
    Returns:
        Coupling coefficient
    """
    # Get segment parameters
    x1 = nec_data.geometry.x[seg1]
    y1 = nec_data.geometry.y[seg1]
    z1 = nec_data.geometry.z[seg1]
    s1 = nec_data.geometry.si[seg1]
    alp1 = nec_data.geometry.alp[seg1]
    bet1 = nec_data.geometry.bet[seg1]
    salp1 = nec_data.geometry.salp[seg1]
    
    x2 = nec_data.geometry.x[seg2]
    y2 = nec_data.geometry.y[seg2]
    z2 = nec_data.geometry.z[seg2]
    s2 = nec_data.geometry.si[seg2]
    alp2 = nec_data.geometry.alp[seg2]
    bet2 = nec_data.geometry.bet[seg2]
    salp2 = nec_data.geometry.salp[seg2]
    
    # Calculate distance between segment centers
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    if r < 1e-10:
        return 0.0
    
    # Calculate coupling using Green's function
    k = 2 * PI / nec_data.geometry.wlam
    green = np.exp(-1j * k * r) / r
    
    # Calculate dot product of segment directions
    dot_product = alp1 * alp2 + bet1 * bet2 + salp1 * salp2
    
    # Calculate coupling coefficient
    coupling = green * dot_product * s1 * s2
    
    return coupling

def calculate_maximum_coupling(y11: complex, y22: complex, y12: complex) -> Tuple[float, complex, complex]:
    """
    Calculate maximum coupling between two segments
    
    Args:
        y11: Self-admittance of segment 1
        y22: Self-admittance of segment 2
        y12: Mutual admittance
        
    Returns:
        Tuple of (coupling_db, load_impedance, input_impedance)
    """
    # Calculate coupling parameter
    yin = y12 * y12
    dbc = abs(yin)
    c = dbc / (2.0 * np.real(y11) * np.real(y22) - np.real(yin))
    
    if c < 0.0 or c > 1.0:
        return -999.99, 0.0, 0.0
    
    # Calculate maximum coupling
    if c < 0.01:
        gmax = 0.5 * (c + 0.25 * c * c * c)
    else:
        gmax = (1.0 - np.sqrt(1.0 - c * c)) / c
    
    # Calculate reflection coefficient
    rho = gmax * np.conj(yin) / dbc
    
    # Calculate load impedance for maximum coupling
    yl = ((1.0 - rho) / (1.0 + rho) + 1.0) * np.real(y22) - y22
    zl = 1.0 / yl
    
    # Calculate input impedance
    yin_total = y11 - yin / (y22 + yl)
    zin = 1.0 / yin_total
    
    # Convert to dB
    coupling_db = db10(gmax)
    
    return coupling_db, zl, zin

def analyze_coupling_network(nec_data: NECData) -> None:
    """
    Analyze coupling for all segments in the network
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\nCOUPLING NETWORK ANALYSIS")
    print("=" * 50)
    
    # Find all segments with the same tag
    tag_groups = {}
    for i in range(nec_data.geometry.n):
        tag = nec_data.geometry.itag[i]
        if tag not in tag_groups:
            tag_groups[tag] = []
        tag_groups[tag].append(i)
    
    # Analyze coupling within each tag group
    for tag, segments in tag_groups.items():
        if len(segments) < 2:
            continue
            
        print(f"\nTag {tag} - {len(segments)} segments:")
        print("-" * 30)
        
        for i in range(len(segments)):
            for j in range(i + 1, len(segments)):
                seg1 = segments[i]
                seg2 = segments[j]
                
                # Calculate coupling
                coupling = calculate_coupling_coefficient(seg1, seg2, nec_data)
                coupling_db = db10(abs(coupling))
                
                print(f"  Segments {seg1+1} - {seg2+1}: {coupling_db:.2f} dB")

def calculate_isolation_matrix(nec_data: NECData) -> np.ndarray:
    """
    Calculate isolation matrix for all segment pairs
    
    Args:
        nec_data: Main NEC data structure
        
    Returns:
        Isolation matrix in dB
    """
    n = nec_data.geometry.n
    isolation = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1, n):
            # Calculate coupling
            coupling = calculate_coupling_coefficient(i, j, nec_data)
            isolation_db = db10(abs(coupling))
            
            isolation[i, j] = isolation_db
            isolation[j, i] = isolation_db
    
    return isolation

def print_isolation_matrix(isolation: np.ndarray, nec_data: NECData) -> None:
    """
    Print isolation matrix
    
    Args:
        isolation: Isolation matrix in dB
        nec_data: Main NEC data structure
    """
    print("\nISOLATION MATRIX (dB)")
    print("=" * 50)
    
    n = min(20, nec_data.geometry.n)  # Print first 20 segments
    
    # Print header
    print("Seg", end="")
    for j in range(n):
        print(f"{j+1:8d}", end="")
    print()
    
    # Print matrix
    for i in range(n):
        print(f"{i+1:3d}", end="")
        for j in range(n):
            if i == j:
                print("     ---", end="")
            else:
                print(f"{isolation[i, j]:8.2f}", end="")
        print()

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

def calculate_coupling_efficiency(nec_data: NECData) -> float:
    """
    Calculate overall coupling efficiency
    
    Args:
        nec_data: Main NEC data structure
        
    Returns:
        Coupling efficiency as percentage
    """
    if nec_data.coupling.ncoup == 0:
        return 0.0
    
    total_coupling = 0.0
    count = 0
    
    for i in range(nec_data.coupling.ncoup):
        if abs(nec_data.coupling.y11a[i]) > 1e-20:
            total_coupling += abs(nec_data.coupling.y11a[i])
            count += 1
    
    if count == 0:
        return 0.0
    
    avg_coupling = total_coupling / count
    efficiency = 100.0 * avg_coupling / abs(nec_data.excitation.vsant[0])
    
    return efficiency

def print_coupling_summary(nec_data: NECData) -> None:
    """
    Print coupling analysis summary
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\nCOUPLING ANALYSIS SUMMARY")
    print("=" * 40)
    print(f"Number of coupling pairs: {nec_data.coupling.ncoup}")
    print(f"Coupling efficiency: {calculate_coupling_efficiency(nec_data):.2f}%")
    
    if nec_data.coupling.ncoup > 0:
        # Find maximum coupling
        max_coupling = -999.99
        max_pair = (0, 0)
        
        for i in range(nec_data.coupling.ncoup):
            coupling_db = db10(abs(nec_data.coupling.y11a[i]))
            if coupling_db > max_coupling:
                max_coupling = coupling_db
                max_pair = (nec_data.coupling.nctag[i], nec_data.coupling.ncseg[i])
        
        print(f"Maximum coupling: {max_coupling:.2f} dB (Tag {max_pair[0]}, Seg {max_pair[1]})")
    
    print("=" * 40) 