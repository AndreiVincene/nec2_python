"""
Matrix Operations Module for Python NEC-2 Implementation

This module contains all matrix-related functions converted from the Fortran code,
including CMSET, FACTRS, SOLVE, and other matrix manipulation subroutines.
"""

import numpy as np
from typing import Tuple, List, Optional
from .constants import *
from .data_structures import *

# ============================================================================
# MAIN MATRIX SETUP (CMSET)
# ============================================================================

def cmset(nec_data: NECData, nrow: int, cm: np.ndarray, rkh: float, iexk: int) -> None:
    """
    Sets up the complex structure matrix in the array CM (CMSET from Fortran)
    
    Args:
        nec_data: Main NEC data structure
        nrow: Number of rows in matrix
        cm: Complex matrix array
        rkh: Integration limit parameter
        iexk: Extended thin wire kernel flag
    """
    mp2 = 2 * nec_data.mp
    npeq = nec_data.np + mp2
    neq = nec_data.n + 2 * nec_data.m
    nop = neq // npeq
    
    if nec_data.icase > 2:
        # Rewind file for out-of-core cases
        pass
    
    nec_data.rkh = rkh
    nec_data.iexk = iexk
    iout = 2 * nec_data.npblk * nrow
    it = nec_data.npblk
    
    # Cycle over matrix blocks
    for ixblk1 in range(1, nec_data.nbloks + 1):
        isv = (ixblk1 - 1) * nec_data.npblk
        if ixblk1 == nec_data.nbloks:
            it = nec_data.nlast
            
        # Initialize matrix block to zero
        cm[:nrow, :it] = 0.0 + 0.0j
        
        i1 = isv + 1
        i2 = isv + it
        in2 = min(i2, nec_data.np)
        im1 = max(i1 - nec_data.np, 1)
        ist = 1
        if i1 <= nec_data.np:
            ist = nec_data.np - i1 + 2
            
        # Wire source loop
        if nec_data.n > 0:
            for j in range(1, nec_data.n + 1):
                trio(j, nec_data)
                
                # Adjust connection indices for matrix structure
                for i in range(nec_data.jsno):
                    ij = nec_data.jco[i]
                    nec_data.jco[i] = ((ij - 1) // nec_data.np) * mp2 + ij
                    
                if i1 <= in2:
                    cmww(j, i1, in2, cm, nrow, cm, nrow, 1, nec_data)
                if im1 <= i2:
                    cmws(j, im1, i2, cm[:, ist-1:], nrow, cm, nrow, 1, nec_data)
                    
                # Matrix elements modified by loading
                if nec_data.nload == 0:
                    continue
                if j > nec_data.np:
                    continue
                ipr = j - isv
                if ipr < 1 or ipr > it:
                    continue
                zaj = nec_data.zarray[j]
                for i in range(nec_data.jsno):
                    jss = nec_data.jco[i]
                    cm[jss-1, ipr-1] -= (nec_data.ax[i] + nec_data.cx[i]) * zaj
                    
        # Matrix elements for patch current sources
        if nec_data.m > 0:
            jm1 = 1 - nec_data.mp
            jm2 = 0
            jst = 1 - mp2
            for i in range(nop):
                jm1 += nec_data.mp
                jm2 += nec_data.mp
                jst += npeq
                if i1 <= in2:
                    cmsw(jm1, jm2, i1, in2, cm[jst-1:jst+npeq-1, :], cm, 0, nrow, 1, nec_data)
                if im1 <= i2:
                    cmss(jm1, jm2, im1, i2, cm[jst-1:jst+npeq-1, ist-1:], nrow, 1, nec_data)
                    
        # Handle symmetry modes
        if nec_data.icase == 1:
            continue
        if nec_data.icase == 3:
            # Write block for out-of-core cases
            blckot(cm, 11, 1, iout, 1, 31)
            continue
            
        # Combine elements for symmetry modes
        for i in range(it):
            for j in range(npeq):
                d = np.zeros(nop, dtype=complex)
                for k in range(nop):
                    ka = j + k * npeq
                    d[k] = cm[ka-1, i]
                deter = np.sum(d)
                cm[j-1, i] = deter
                for k in range(1, nop):
                    ka = j + k * npeq
                    deter = d[0]
                    for kk in range(1, nop):
                        deter += d[kk] * nec_data.ssx[k, kk]
                    cm[ka-1, i] = deter
                    
    if nec_data.icase > 2:
        # Rewind file
        pass

# ============================================================================
# MATRIX FACTORIZATION (FACTRS)
# ============================================================================

def factrs(npeq: int, neq: int, cm: np.ndarray, ip: np.ndarray, ix: np.ndarray, 
           i11: int, i12: int, i13: int, i14: int) -> None:
    """
    Factors the complex structure matrix (FACTRS from Fortran)
    
    Args:
        npeq: Number of equations in primary matrix
        neq: Total number of equations
        cm: Complex matrix array
        ip: Pivot array
        ix: Index array
        i11, i12, i13, i14: File unit numbers
    """
    # Implementation of LU decomposition with pivoting
    # This is a complex matrix factorization routine
    
    for i in range(1, neq + 1):
        ip[i-1] = i
        
    for k in range(1, neq):
        # Find pivot element
        pivot = abs(cm[k-1, k-1])
        ipivot = k
        for i in range(k + 1, neq + 1):
            if abs(cm[i-1, k-1]) > pivot:
                pivot = abs(cm[i-1, k-1])
                ipivot = i
                
        if ipivot != k:
            # Swap rows
            ip[k-1], ip[ipivot-1] = ip[ipivot-1], ip[k-1]
            for j in range(k, neq + 1):
                cm[k-1, j-1], cm[ipivot-1, j-1] = cm[ipivot-1, j-1], cm[k-1, j-1]
                
        # Eliminate column
        for i in range(k + 1, neq + 1):
            factor = cm[i-1, k-1] / cm[k-1, k-1]
            cm[i-1, k-1] = factor
            for j in range(k + 1, neq + 1):
                cm[i-1, j-1] -= factor * cm[k-1, j-1]

# ============================================================================
# MATRIX SOLVER (SOLVE)
# ============================================================================

def solve(npeq: int, neq: int, cm: np.ndarray, ip: np.ndarray, ix: np.ndarray,
          b: np.ndarray, i11: int, i12: int, i13: int, i14: int) -> None:
    """
    Solves the complex matrix equation (SOLVE from Fortran)
    
    Args:
        npeq: Number of equations in primary matrix
        neq: Total number of equations
        cm: Complex matrix array (factored)
        ip: Pivot array
        ix: Index array
        b: Right-hand side vector
        i11, i12, i13, i14: File unit numbers
    """
    # Forward substitution
    for i in range(1, neq + 1):
        ipivot = ip[i-1]
        if ipivot != i:
            b[i-1], b[ipivot-1] = b[ipivot-1], b[i-1]
            
    for i in range(2, neq + 1):
        for j in range(1, i):
            b[i-1] -= cm[i-1, j-1] * b[j-1]
            
    # Backward substitution
    b[neq-1] /= cm[neq-1, neq-1]
    for i in range(neq - 1, 0, -1):
        for j in range(i + 1, neq + 1):
            b[i-1] -= cm[i-1, j-1] * b[j-1]
        b[i-1] /= cm[i-1, i-1]

# ============================================================================
# WIRE-WIRE MATRIX ELEMENTS (CMWW)
# ============================================================================

def cmww(j: int, i1: int, i2: int, cm: np.ndarray, nr: int, cw: np.ndarray, 
         nw: int, itrp: int, nec_data: NECData) -> None:
    """
    Computes matrix elements for wire-wire interactions (CMWW from Fortran)
    
    Args:
        j: Source segment number
        i1, i2: Observation segment range
        cm: Matrix array
        nr: Number of rows
        cw: Additional matrix array
        nw: Number of columns in cw
        itrp: Transpose flag
        nec_data: Main NEC data structure
    """
    # Set source segment parameters
    s = nec_data.si[j-1]
    b = nec_data.bi[j-1]
    xj = nec_data.x[j-1]
    yj = nec_data.y[j-1]
    zj = nec_data.z[j-1]
    cabj = nec_data.alp[j-1]
    sabj = nec_data.bet[j-1]
    salpj = nec_data.salp[j-1]
    
    # Determine extended thin wire approximation flags
    if nec_data.iexk == 0:
        ind1 = 0
        ind2 = 0
    else:
        # Check connection conditions for extended thin wire approximation
        ipr = nec_data.icon1[j-1]
        if ipr > 10000:
            ind1 = 0
        elif ipr < 0:
            ipr = -ipr
            if -nec_data.icon1[ipr-1] == j:
                ind1 = 0
            else:
                ind1 = 2
        elif ipr == j:
            if cabj**2 + sabj**2 > 1e-8:
                ind1 = 2
            else:
                ind1 = 0
        elif nec_data.icon2[ipr-1] == j:
            xi = abs(cabj * nec_data.alp[ipr-1] + sabj * nec_data.bet[ipr-1] + 
                    salpj * nec_data.salp[ipr-1])
            if xi < 0.999999:
                ind1 = 2
            elif abs(nec_data.bi[ipr-1] / b - 1.0) > 1e-6:
                ind1 = 2
            else:
                ind1 = 0
        else:
            ind1 = 2
            
        # Similar check for end 2
        ipr = nec_data.icon2[j-1]
        if ipr > 10000:
            ind2 = 0
        elif ipr < 0:
            ipr = -ipr
            if -nec_data.icon2[ipr-1] == j:
                ind2 = 0
            else:
                ind2 = 2
        elif ipr == j:
            if cabj**2 + sabj**2 > 1e-8:
                ind2 = 2
            else:
                ind2 = 0
        elif nec_data.icon1[ipr-1] == j:
            xi = abs(cabj * nec_data.alp[ipr-1] + sabj * nec_data.bet[ipr-1] + 
                    salpj * nec_data.salp[ipr-1])
            if xi < 0.999999:
                ind2 = 2
            elif abs(nec_data.bi[ipr-1] / b - 1.0) > 1e-6:
                ind2 = 2
            else:
                ind2 = 0
        else:
            ind2 = 2
    
    # Observation loop
    ipr = 0
    for i in range(i1, i2 + 1):
        ipr += 1
        ij = i - j
        xi = nec_data.x[i-1]
        yi = nec_data.y[i-1]
        zi = nec_data.z[i-1]
        ai = nec_data.bi[i-1]
        cabi = nec_data.alp[i-1]
        sabi = nec_data.bet[i-1]
        salpi = nec_data.salp[i-1]
        
        # Calculate field components
        efield(xi, yi, zi, ai, ij, nec_data)
        etk = (nec_data.exk * cabi + nec_data.eyk * sabi + nec_data.ezk * salpi)
        ets = (nec_data.exs * cabi + nec_data.eys * sabi + nec_data.ezs * salpi)
        etc = (nec_data.exc * cabi + nec_data.eyc * sabi + nec_data.ezc * salpi)
        
        # Fill matrix elements
        if itrp == 0:
            # Normal fill
            for ij in range(nec_data.jsno):
                jx = nec_data.jco[ij]
                cm[ipr-1, jx-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                   etc * nec_data.cx[ij])
        elif itrp == 2:
            # Transposed fill for C(WW) - test for elements for D(WW)PRIME
            for ij in range(nec_data.jsno):
                jx = nec_data.jco[ij]
                if jx <= nr:
                    cm[jx-1, ipr-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                       etc * nec_data.cx[ij])
                else:
                    jx = jx - nr
                    cw[jx-1, ipr-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                       etc * nec_data.cx[ij])
        else:
            # Transposed fill
            for ij in range(nec_data.jsno):
                jx = nec_data.jco[ij]
                cm[jx-1, ipr-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                   etc * nec_data.cx[ij])

# ============================================================================
# WIRE-SURFACE MATRIX ELEMENTS (CMWS)
# ============================================================================

def cmws(j: int, i1: int, i2: int, cm: np.ndarray, nr: int, cw: np.ndarray, 
         nw: int, itrp: int, nec_data: NECData) -> None:
    """
    Computes matrix elements for wire-surface interactions (CMWS from Fortran)
    
    Args:
        j: Source segment number
        i1, i2: Observation segment range
        cm: Matrix array
        nr: Number of rows
        cw: Additional matrix array
        nw: Number of columns in cw
        itrp: Transpose flag
        nec_data: Main NEC data structure
    """
    # Set source segment parameters
    s = nec_data.si[j-1]
    b = nec_data.bi[j-1]
    xj = nec_data.x[j-1]
    yj = nec_data.y[j-1]
    zj = nec_data.z[j-1]
    cabj = nec_data.alp[j-1]
    sabj = nec_data.bet[j-1]
    salpj = nec_data.salp[j-1]
    
    # Observation loop
    ipr = 0
    for i in range(i1, i2 + 1):
        ipr += 1
        ipatch = (i + 1) // 2
        ik = i - (i // 2) * 2
        
        if ik == 0 and ipr != 1:
            js = nec_data.ld + 1 - ipatch
            xi = nec_data.x[js-1]
            yi = nec_data.y[js-1]
            zi = nec_data.z[js-1]
            hsfield(xi, yi, zi, 0.0, nec_data)
            if ik == 0:
                tx = nec_data.t2x[js-1]
                ty = nec_data.t2y[js-1]
                tz = nec_data.t2z[js-1]
            else:
                tx = nec_data.t1x[js-1]
                ty = nec_data.t1y[js-1]
                tz = nec_data.t1z[js-1]
        else:
            js = nec_data.ld + 1 - ipatch
            xi = nec_data.x[js-1]
            yi = nec_data.y[js-1]
            zi = nec_data.z[js-1]
            hsfield(xi, yi, zi, 0.0, nec_data)
            tx = nec_data.t1x[js-1]
            ty = nec_data.t1y[js-1]
            tz = nec_data.t1z[js-1]
            
        etk = -(nec_data.exk * tx + nec_data.eyk * ty + nec_data.ezk * tz) * nec_data.salp[js-1]
        ets = -(nec_data.exs * tx + nec_data.eys * ty + nec_data.ezs * tz) * nec_data.salp[js-1]
        etc = -(nec_data.exc * tx + nec_data.eyc * ty + nec_data.ezc * tz) * nec_data.salp[js-1]
        
        # Fill matrix elements
        if itrp == 0:
            # Normal fill
            for ij in range(nec_data.jsno):
                jx = nec_data.jco[ij]
                cm[ipr-1, jx-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                   etc * nec_data.cx[ij])
        elif itrp == 2:
            # Transposed fill - C(WS) and D(WS)PRIME (=CW)
            for ij in range(nec_data.jsno):
                jx = nec_data.jco[ij]
                if jx <= nr:
                    cm[jx-1, ipr-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                       etc * nec_data.cx[ij])
                else:
                    jx = jx - nr
                    cw[jx-1, ipr-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                       etc * nec_data.cx[ij])
        else:
            # Transposed fill
            for ij in range(nec_data.jsno):
                jx = nec_data.jco[ij]
                cm[jx-1, ipr-1] += (etk * nec_data.ax[ij] + ets * nec_data.bx[ij] + 
                                   etc * nec_data.cx[ij])

# ============================================================================
# SURFACE-SURFACE MATRIX ELEMENTS (CMSS)
# ============================================================================

def cmss(j1: int, j2: int, im1: int, im2: int, cm: np.ndarray, nrow: int, 
         itrp: int, nec_data: NECData) -> None:
    """
    Computes matrix elements for surface-surface interactions (CMSS from Fortran)
    
    Args:
        j1, j2: Source patch range
        im1, im2: Observation patch range
        cm: Matrix array
        nrow: Number of rows
        itrp: Transpose flag
        nec_data: Main NEC data structure
    """
    ldp = nec_data.ld + 1
    i1 = (im1 + 1) // 2
    i2 = (im2 + 1) // 2
    icomp = i1 * 2 - 3
    ii1 = -1
    if icomp + 2 < im1:
        ii1 = -2
        
    # Loop over observation patches
    for i in range(i1, i2 + 1):
        il = ldp - i
        icomp += 2
        ii1 += 2
        ii2 = ii1 + 1
        
        t1xi = nec_data.t1x[il-1] * nec_data.salp[il-1]
        t1yi = nec_data.t1y[il-1] * nec_data.salp[il-1]
        t1zi = nec_data.t1z[il-1] * nec_data.salp[il-1]
        t2xi = nec_data.t2x[il-1] * nec_data.salp[il-1]
        t2yi = nec_data.t2y[il-1] * nec_data.salp[il-1]
        t2zi = nec_data.t2z[il-1] * nec_data.salp[il-1]
        xi = nec_data.x[il-1]
        yi = nec_data.y[il-1]
        zi = nec_data.z[il-1]
        jj1 = -1
        
        # Loop over source patches
        for j in range(j1, j2 + 1):
            jl = ldp - j
            jj1 += 2
            jj2 = jj1 + 1
            s = nec_data.bi[jl-1]
            xj = nec_data.x[jl-1]
            yj = nec_data.y[jl-1]
            zj = nec_data.z[jl-1]
            t1xj = nec_data.t1x[jl-1]
            t1yj = nec_data.t1y[jl-1]
            t1zj = nec_data.t1z[jl-1]
            t2xj = nec_data.t2x[jl-1]
            t2yj = nec_data.t2y[jl-1]
            t2zj = nec_data.t2z[jl-1]
            
            # Calculate field components
            hintg(xi, yi, zi, nec_data)
            g11 = -(t2xi * nec_data.exk + t2yi * nec_data.eyk + t2zi * nec_data.ezk)
            g12 = -(t2xi * nec_data.exs + t2yi * nec_data.eys + t2zi * nec_data.ezs)
            g21 = -(t1xi * nec_data.exk + t1yi * nec_data.eyk + t1zi * nec_data.ezk)
            g22 = -(t1xi * nec_data.exs + t1yi * nec_data.eys + t1zi * nec_data.ezs)
            
            if i != j:
                g11 -= 0.5
                g22 += 0.5
                
            if itrp == 0:
                # Normal fill
                if icomp >= im1:
                    cm[ii1-1, jj1-1] = g11
                    cm[ii1-1, jj2-1] = g12
                if icomp < im2:
                    cm[ii2-1, jj1-1] = g21
                    cm[ii2-1, jj2-1] = g22
            else:
                # Transposed fill
                if icomp >= im1:
                    cm[jj1-1, ii1-1] = g11
                    cm[jj2-1, ii1-1] = g12
                if icomp < im2:
                    cm[jj1-1, ii2-1] = g21
                    cm[jj2-1, ii2-1] = g22

# ============================================================================
# SURFACE-WIRE MATRIX ELEMENTS (CMSW)
# ============================================================================

def cmsw(j1: int, j2: int, i1: int, i2: int, cm: np.ndarray, cw: np.ndarray, 
         ncw: int, nrow: int, itrp: int, nec_data: NECData) -> None:
    """
    Computes matrix elements for E along wires due to patch current (CMSW from Fortran)
    
    Args:
        j1, j2: Source patch range
        i1, i2: Observation segment range
        cm: Matrix array
        cw: Additional matrix array
        ncw: Number of columns in cw
        nrow: Number of rows
        itrp: Transpose flag
        nec_data: Main NEC data structure
    """
    # Implementation of surface-wire matrix element calculation
    # This is a complex routine that handles patch-to-wire interactions
    
    ldp = nec_data.ld + 1
    neqs = nec_data.n - nec_data.n1 + 2 * (nec_data.m - nec_data.m1)
    
    if itrp < 0:
        # Handle special case for old segments connecting to patches
        return
        
    k = 0
    icgo = 1
    
    # Observation loop
    for i in range(i1, i2 + 1):
        k += 1
        xi = nec_data.x[i-1]
        yi = nec_data.y[i-1]
        zi = nec_data.z[i-1]
        cabi = nec_data.alp[i-1]
        sabi = nec_data.bet[i-1]
        salpi = nec_data.salp[i-1]
        
        ipch = 0
        if nec_data.icon1[i-1] < 10000:
            pass
        else:
            ipch = nec_data.icon1[i-1] - 10000
            fsign = -1.0
            
        if nec_data.icon2[i-1] < 10000:
            pass
        else:
            ipch = nec_data.icon2[i-1] - 10000
            fsign = 1.0
            
        jl = 0
        
        # Source loop
        for j in range(j1, j2 + 1):
            js = ldp - j
            jl += 2
            t1xj = nec_data.t1x[js-1]
            t1yj = nec_data.t1y[js-1]
            t1zj = nec_data.t1z[js-1]
            t2xj = nec_data.t2x[js-1]
            t2yj = nec_data.t2y[js-1]
            t2zj = nec_data.t2z[js-1]
            xj = nec_data.x[js-1]
            yj = nec_data.y[js-1]
            zj = nec_data.z[js-1]
            s = nec_data.bi[js-1]
            
            # Ground loop
            for ip in range(1, nec_data.ksymp + 1):
                nec_data.ipgnd = ip
                
                if ipch != j and icgo == 1:
                    # Calculate field using UNERE
                    unere(xi, yi, zi, nec_data)
                    if itrp == 0:
                        # Normal fill
                        cm[k-1, jl-2] += (nec_data.exk * cabi + nec_data.eyk * sabi + 
                                         nec_data.ezk * salpi)
                        cm[k-1, jl-1] += (nec_data.exs * cabi + nec_data.eys * sabi + 
                                         nec_data.ezs * salpi)
                    else:
                        # Transposed fill
                        cm[jl-2, k-1] += (nec_data.exk * cabi + nec_data.eyk * sabi + 
                                         nec_data.ezk * salpi)
                        cm[jl-1, k-1] += (nec_data.exs * cabi + nec_data.eys * sabi + 
                                         nec_data.ezs * salpi)
                else:
                    # Handle patch connection case
                    if ip == 2:
                        continue
                    if icgo > 1:
                        pass
                    else:
                        # Calculate patch current integration
                        pcint(xi, yi, zi, cabi, sabi, salpi, nec_data)
                        py = np.pi * nec_data.si[i-1] * fsign
                        px = np.sin(py)
                        py = np.cos(py)
                        exc = nec_data.emel[8] * fsign  # Index 8 for 9th element
                        trio(i, nec_data)
                        
                        if i > nec_data.n1:
                            il = i - ncw
                            if i <= nec_data.np:
                                il = ((il - 1) // nec_data.np) * 2 * nec_data.mp + il
                        else:
                            il = neqs + nec_data.iconx[i-1]
                            
                        if itrp == 0:
                            cw[k-1, il-1] += exc * (nec_data.ax[nec_data.jsno-1] + 
                                                   nec_data.bx[nec_data.jsno-1] * px + 
                                                   nec_data.cx[nec_data.jsno-1] * py)
                        else:
                            cw[il-1, k-1] += exc * (nec_data.ax[nec_data.jsno-1] + 
                                                   nec_data.bx[nec_data.jsno-1] * px + 
                                                   nec_data.cx[nec_data.jsno-1] * py)
                            
                    if itrp == 0:
                        cm[k-1, jl-2] = nec_data.emel[icgo-1]
                        cm[k-1, jl-1] = nec_data.emel[icgo+3]  # icgo+4-1
                    else:
                        cm[jl-2, k-1] = nec_data.emel[icgo-1]
                        cm[jl-1, k-1] = nec_data.emel[icgo+3]
                        
                    icgo += 1
                    if icgo == 5:
                        icgo = 1

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def trio(j: int, nec_data: NECData) -> None:
    """
    Sets up triangle functions for segment j (TRIO from Fortran)
    
    Args:
        j: Segment number
        nec_data: Main NEC data structure
    """
    # Implementation of triangle function setup
    # This calculates the basis functions for segment j
    
    # Get segment parameters
    s = nec_data.si[j-1]
    b = nec_data.bi[j-1]
    xj = nec_data.x[j-1]
    yj = nec_data.y[j-1]
    zj = nec_data.z[j-1]
    cabj = nec_data.alp[j-1]
    sabj = nec_data.bet[j-1]
    salpj = nec_data.salp[j-1]
    
    # Calculate connection data
    jsno = 0
    jco = []
    ax = []
    bx = []
    cx = []
    
    # Check end 1 connection
    icon1 = nec_data.icon1[j-1]
    if icon1 != 0:
        jsno += 1
        if icon1 < 0:
            jco.append(-icon1)
            ax.append(1.0)
            bx.append(0.0)
            cx.append(0.0)
        else:
            jco.append(icon1)
            ax.append(1.0)
            bx.append(0.0)
            cx.append(0.0)
            
    # Check end 2 connection
    icon2 = nec_data.icon2[j-1]
    if icon2 != 0:
        jsno += 1
        if icon2 < 0:
            jco.append(-icon2)
            ax.append(0.0)
            bx.append(0.0)
            cx.append(1.0)
        else:
            jco.append(icon2)
            ax.append(0.0)
            bx.append(0.0)
            cx.append(1.0)
            
    # Add self term
    jsno += 1
    jco.append(j)
    ax.append(0.0)
    bx.append(1.0)
    cx.append(0.0)
    
    # Store results
    nec_data.jsno = jsno
    nec_data.jco[:jsno] = jco
    nec_data.ax[:jsno] = ax
    nec_data.bx[:jsno] = bx
    nec_data.cx[:jsno] = cx

def efield(xi: float, yi: float, zi: float, ai: float, ij: int, nec_data: NECData) -> None:
    """
    Compute near E fields of a segment (EFLD from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        ai: Observation segment radius
        ij: Source-observation segment difference
        nec_data: Main NEC data structure
    """
    # This is a placeholder for the complex field calculation
    # The actual implementation would calculate EXK, EYK, EZK, etc.
    
    # Set default values (these would be calculated in the real implementation)
    nec_data.exk = 0.0 + 0.0j
    nec_data.eyk = 0.0 + 0.0j
    nec_data.ezk = 0.0 + 0.0j
    nec_data.exs = 0.0 + 0.0j
    nec_data.eys = 0.0 + 0.0j
    nec_data.ezs = 0.0 + 0.0j
    nec_data.exc = 0.0 + 0.0j
    nec_data.eyc = 0.0 + 0.0j
    nec_data.ezc = 0.0 + 0.0j

def hsfield(xi: float, yi: float, zi: float, ai: float, nec_data: NECData) -> None:
    """
    Compute H field components (HSFLD from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        ai: Observation segment radius
        nec_data: Main NEC data structure
    """
    # Placeholder for H field calculation
    pass

def hintg(xi: float, yi: float, zi: float, nec_data: NECData) -> None:
    """
    Compute H field integration (HINTG from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        nec_data: Main NEC data structure
    """
    # Placeholder for H field integration
    pass

def unere(xi: float, yi: float, zi: float, nec_data: NECData) -> None:
    """
    Compute E field using UNERE (UNERE from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        nec_data: Main NEC data structure
    """
    # Placeholder for UNERE field calculation
    pass

def pcint(xi: float, yi: float, zi: float, cabi: float, sabi: float, 
          salpi: float, nec_data: NECData) -> None:
    """
    Compute patch current integration (PCINT from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        cabi, sabi, salpi: Direction cosines
        nec_data: Main NEC data structure
    """
    # Placeholder for patch current integration
    # Initialize emel array
    nec_data.emel = np.zeros(9, dtype=complex)

def blckot(ar: np.ndarray, nunit: int, ix1: int, ix2: int, nblks: int, neof: int) -> None:
    """
    Controls reading and writing of matrix blocks on files (BLCKOT from Fortran)
    
    Args:
        ar: Array to write
        nunit: File unit number
        ix1, ix2: Index range
        nblks: Number of blocks
        neof: EOF flag
    """
    # Placeholder for block output
    pass 