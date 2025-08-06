"""
Ground Effects Module for Python NEC-2 Implementation

This module contains all ground-related functions converted from the Fortran code,
including Sommerfeld calculations, ground reflection coefficients, and other ground effect subroutines.
"""

import numpy as np
from typing import Tuple, List, Optional, Complex
from .constants import *
from .data_structures import *

# ============================================================================
# SOMMERFELD GROUND CALCULATIONS (SOM2D)
# ============================================================================

def som2d(fmhz: float, epr: float, sig: float) -> None:
    """
    Compute Sommerfeld ground integral tables (SOM2D from Fortran)
    
    Args:
        fmhz: Frequency in MHz
        epr: Relative dielectric constant
        sig: Conductivity in mhos/meter
    """
    # Initialize complex constants
    ck2 = 2 * PI
    ck2sq = ck2 * ck2
    
    # Calculate wavelength and complex dielectric constant
    wlam = CVEL / fmhz
    if sig < 0:
        epscf = complex(epr, sig)
    else:
        epscf = complex(epr, -sig * wlam * 59.96)
    
    # Calculate wave numbers
    ck1sq = ck2sq * np.conj(epscf)
    ck1 = np.sqrt(ck1sq)
    ck1r = np.real(ck1)
    
    # Calculate magnitude parameters
    tkmag = 100.0 * abs(ck1)
    tsmag = 100.0 * ck1 * np.conj(ck1)
    
    # Calculate Sommerfeld parameters
    cksm = ck2sq / (ck1sq + ck2sq)
    ct1 = 0.5 * (ck1sq - ck2sq)
    erv = ck1sq * ck1sq
    ezv = ck2sq * ck2sq
    ct2 = 0.125 * (erv - ezv)
    erv = erv * ck1sq
    ezv = ezv * ck2sq
    ct3 = 0.0625 * (erv - ezv)
    
    # Initialize grid arrays
    ar1 = np.zeros((11, 10, 4), dtype=complex)
    ar2 = np.zeros((17, 5, 4), dtype=complex)
    ar3 = np.zeros((9, 8, 4), dtype=complex)
    
    # Grid parameters
    nxa = [11, 17, 9]
    nya = [10, 5, 8]
    xsa = [0.0, 0.2, 0.2]
    ysa = [0.0, 0.0, 0.3490658504]
    dxa = [0.02, 0.05, 0.1]
    dya = [0.1745329252, 0.0872664626, 0.1745329252]
    
    # Loop over 3 grid regions
    for k in range(3):
        nr = nxa[k]
        nth = nya[k]
        dr = dxa[k]
        dth = dya[k]
        r = xsa[k] - dr
        irs = 1
        if k == 0:
            r = xsa[k]
            irs = 2
            
        # Loop over R
        for ir in range(irs, nr + 1):
            r += dr
            thet = ysa[k] - dth
            
            # Loop over theta
            for ith in range(nth):
                thet += dth
                rho = r * np.cos(thet)
                zph = r * np.sin(thet)
                
                if rho < 1e-7:
                    rho = 1e-8
                if zph < 1e-7:
                    zph = 0.0
                    
                # Calculate Sommerfeld integrals
                erv_val, ezv_val, erh_val, eph_val = evlua(rho, zph, ck1, ck2, ck1sq, ck2sq)
                
                rk = ck2 * r
                con = -4.77147j * r / complex(np.cos(rk), -np.sin(rk))
                
                # Store results in appropriate grid
                if k == 0:
                    ar1[ir-1, ith, 0] = erv_val * con
                    ar1[ir-1, ith, 1] = ezv_val * con
                    ar1[ir-1, ith, 2] = erh_val * con
                    ar1[ir-1, ith, 3] = eph_val * con
                elif k == 1:
                    ar2[ir-1, ith, 0] = erv_val * con
                    ar2[ir-1, ith, 1] = ezv_val * con
                    ar2[ir-1, ith, 2] = erh_val * con
                    ar2[ir-1, ith, 3] = eph_val * con
                else:
                    ar3[ir-1, ith, 0] = erv_val * con
                    ar3[ir-1, ith, 1] = ezv_val * con
                    ar3[ir-1, ith, 2] = erh_val * con
                    ar3[ir-1, ith, 3] = eph_val * con
    
    # Fill grid 1 for R equal to zero
    cl2 = -188.370j * (epscf - 1.0) / (epscf + 1.0)
    cl1 = cl2 / (epscf + 1.0)
    ezv_val = epscf * cl1
    
    thet = -dya[0]
    nth = nya[0]
    for ith in range(nth):
        thet += dya[0]
        if ith == nth - 1:
            erv_val = 0.0
            erh_val = cl2 - 0.5 * cl1
            eph_val = -erh_val
        else:
            tfac2 = np.cos(thet)
            tfac1 = (1.0 - np.sin(thet)) / tfac2
            tfac2 = tfac1 / tfac2
            erv_val = epscf * cl1 * tfac1
            erh_val = cl1 * (tfac2 - 1.0) + cl2
            eph_val = cl1 * tfac2 - cl2
            
        ar1[0, ith, 0] = erv_val
        ar1[0, ith, 1] = ezv_val
        ar1[0, ith, 2] = erh_val
        ar1[0, ith, 3] = eph_val

def evlua(rho: float, zph: float, ck1: complex, ck2: float, 
          ck1sq: complex, ck2sq: float) -> Tuple[complex, complex, complex, complex]:
    """
    Evaluate Sommerfeld integrals (EVLUA from Fortran)
    
    Args:
        rho: Radial distance
        zph: Z coordinate
        ck1: Complex wave number 1
        ck2: Wave number 2
        ck1sq: Complex wave number 1 squared
        ck2sq: Wave number 2 squared
        
    Returns:
        erv, ezv, erh, eph: Sommerfeld integral values
    """
    # Set up integration parameters
    del_val = zph
    if rho > del_val:
        del_val = rho
    if zph < 2.0 * rho:
        # Use Hankel function form
        jh = 1
        cp1 = complex(0.0, 0.4 * ck2)
        cp2 = complex(0.6 * ck2, -0.2 * ck2)
        cp3 = complex(1.02 * ck2, -0.2 * ck2)
        
        # Initialize integration
        a = cp1
        b = cp2
        sum_vals = rom1_hankel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
        
        a = cp2
        b = cp3
        ans_vals = rom1_hankel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
        
        for i in range(6):
            sum_vals[i] = -(sum_vals[i] + ans_vals[i])
            
        # Path from imaginary axis to -infinity
        slope = 1000.0
        if zph > 0.001 * rho:
            slope = rho / zph
            
        del_val = PI / del_val
        delta = complex(-1.0, slope) * del_val / np.sqrt(1.0 + slope * slope)
        delta2 = -np.conj(delta)
        
        # Integrate using GSHANK
        ans_vals = gshank(cp1, delta, sum_vals, 6, rho, zph, ck1, ck2, ck1sq, ck2sq)
        
    else:
        # Use Bessel function form
        jh = 0
        a = complex(0.0, 0.0)
        del_val = 1.0 / del_val
        
        if del_val <= abs(ck1):
            b = complex(del_val, -del_val)
            ans_vals = rom1_bessel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
        else:
            b = complex(0.1 * abs(ck1), -0.1 * abs(ck1))
            sum_vals = rom1_bessel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
            a = b
            b = complex(del_val, -del_val)
            ans_vals = rom1_bessel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
            for i in range(6):
                ans_vals[i] = sum_vals[i] + ans_vals[i]
                
        delta = PI * del_val
        ans_vals = gshank(b, delta, ans_vals, 6, rho, zph, ck1, ck2, ck1sq, ck2sq)
    
    # Final calculations
    ans_vals[5] = ans_vals[5] * ck1
    
    # Conjugate since NEC uses exp(+jwt)
    erv = np.conj(ck1sq * ans_vals[2])
    ezv = np.conj(ck1sq * (ans_vals[1] + ck2sq * ans_vals[4]))
    erh = np.conj(ck2sq * (ans_vals[0] + ans_vals[5]))
    eph = -np.conj(ck2sq * (ans_vals[3] + ans_vals[5]))
    
    return erv, ezv, erh, eph

def rom1_bessel(a: complex, b: complex, rho: float, zph: float,
                ck1: complex, ck2: float, ck1sq: complex, ck2sq: float) -> List[complex]:
    """
    Romberg integration for Bessel function form (ROM1 from Fortran)
    
    Args:
        a, b: Integration limits
        rho, zph: Coordinates
        ck1, ck2: Wave numbers
        ck1sq, ck2sq: Wave numbers squared
        
    Returns:
        List of 6 integral values
    """
    # Initialize result array
    sum_vals = [complex(0.0, 0.0)] * 6
    
    # Integration parameters
    nm = 131072
    nts = 4
    rx = 1e-4
    
    # Perform Romberg integration
    # This is a simplified implementation
    # The full implementation would use adaptive quadrature
    
    # For now, use a simple trapezoidal rule
    n_points = 100
    dt = (b - a) / n_points
    
    for i in range(n_points + 1):
        t = a + i * dt
        if i == 0 or i == n_points:
            weight = 0.5
        else:
            weight = 1.0
            
        # Calculate integrand
        integrand = saoa_bessel(t, rho, zph, ck1, ck2, ck1sq, ck2sq)
        
        for j in range(6):
            sum_vals[j] += weight * dt * integrand[j]
    
    return sum_vals

def rom1_hankel(a: complex, b: complex, rho: float, zph: float,
                ck1: complex, ck2: float, ck1sq: complex, ck2sq: float) -> List[complex]:
    """
    Romberg integration for Hankel function form (ROM1 from Fortran)
    
    Args:
        a, b: Integration limits
        rho, zph: Coordinates
        ck1, ck2: Wave numbers
        ck1sq, ck2sq: Wave numbers squared
        
    Returns:
        List of 6 integral values
    """
    # Initialize result array
    sum_vals = [complex(0.0, 0.0)] * 6
    
    # Integration parameters
    nm = 131072
    nts = 4
    rx = 1e-4
    
    # Perform Romberg integration
    # This is a simplified implementation
    # The full implementation would use adaptive quadrature
    
    # For now, use a simple trapezoidal rule
    n_points = 100
    dt = (b - a) / n_points
    
    for i in range(n_points + 1):
        t = a + i * dt
        if i == 0 or i == n_points:
            weight = 0.5
        else:
            weight = 1.0
            
        # Calculate integrand
        integrand = saoa_hankel(t, rho, zph, ck1, ck2, ck1sq, ck2sq)
        
        for j in range(6):
            sum_vals[j] += weight * dt * integrand[j]
    
    return sum_vals

def saoa_bessel(t: complex, rho: float, zph: float, ck1: complex, ck2: float,
                ck1sq: complex, ck2sq: float) -> List[complex]:
    """
    Calculate integrand for Bessel function form (SAOA from Fortran)
    
    Args:
        t: Integration parameter
        rho, zph: Coordinates
        ck1, ck2: Wave numbers
        ck1sq, ck2sq: Wave numbers squared
        
    Returns:
        List of 6 integrand values
    """
    # Calculate lambda from parameter t
    xlam = t  # Simplified - full implementation would use LAMBDA function
    
    # Calculate Bessel functions
    b0, b0p = bessel(xlam * rho)
    b0 = 2.0 * b0
    b0p = 2.0 * b0p
    
    # Calculate gamma functions
    cgam1 = np.sqrt(xlam * xlam - ck1sq)
    cgam2 = np.sqrt(xlam * xlam - ck2sq)
    
    if np.real(cgam1) == 0.0:
        cgam1 = complex(0.0, -abs(np.imag(cgam1)))
    if np.real(cgam2) == 0.0:
        cgam2 = complex(0.0, -abs(np.imag(cgam2)))
    
    # Calculate integrand components
    xlr = xlam * np.conj(xlam)
    if xlr < abs(ck1 * np.conj(ck1)):
        dgam = cgam2 - cgam1
    else:
        if np.imag(xlam) < 0.0:
            sign = 1.0
        else:
            xlr = np.real(xlam)
            if xlr < ck2:
                sign = -1.0
            elif xlr > np.real(ck1):
                sign = 1.0
            else:
                dgam = cgam2 - cgam1
                sign = 0.0
                
        if sign != 0.0:
            dgam = 1.0 / (xlam * xlam)
            dgam = sign * ((ct3 * dgam + ct2) * dgam + ct1) / xlam
    
    # Calculate denominators
    den2 = cksm * dgam / (cgam2 * (ck1sq * cgam2 + ck2sq * cgam1))
    den1 = 1.0 / (cgam1 + cgam2) - cksm / cgam2
    
    # Calculate exponential factor
    com = np.exp(-cgam2 * zph)
    
    # Calculate results
    ans = [complex(0.0, 0.0)] * 6
    ans[5] = com * b0 * den1 / ck1
    com = com * den2
    
    if rho == 0.0:
        ans[0] = -com * xlam * xlam * 0.5
        ans[3] = ans[0]
    else:
        b0p = b0p / rho
        ans[0] = -com * xlam * (b0p + b0 * xlam)
        ans[3] = com * xlam * b0p
        
    ans[1] = com * cgam2 * cgam2 * b0
    ans[2] = -ans[3] * cgam2 * rho
    ans[4] = com * b0
    
    return ans

def saoa_hankel(t: complex, rho: float, zph: float, ck1: complex, ck2: float,
                ck1sq: complex, ck2sq: float) -> List[complex]:
    """
    Calculate integrand for Hankel function form (SAOA from Fortran)
    
    Args:
        t: Integration parameter
        rho, zph: Coordinates
        ck1, ck2: Wave numbers
        ck1sq, ck2sq: Wave numbers squared
        
    Returns:
        List of 6 integrand values
    """
    # Calculate lambda from parameter t
    xlam = t  # Simplified - full implementation would use LAMBDA function
    
    # Calculate Hankel functions
    b0, b0p = hankel(xlam * rho)
    
    # Calculate gamma functions
    com = xlam - ck1
    cgam1 = np.sqrt(xlam + ck1) * np.sqrt(com)
    if np.real(com) < 0.0 and np.imag(com) >= 0.0:
        cgam1 = -cgam1
        
    com = xlam - ck2
    cgam2 = np.sqrt(xlam + ck2) * np.sqrt(com)
    if np.real(com) < 0.0 and np.imag(com) >= 0.0:
        cgam2 = -cgam2
    
    # Calculate integrand components
    xlr = xlam * np.conj(xlam)
    if xlr < abs(ck1 * np.conj(ck1)):
        dgam = cgam2 - cgam1
    else:
        if np.imag(xlam) < 0.0:
            sign = 1.0
        else:
            xlr = np.real(xlam)
            if xlr < ck2:
                sign = -1.0
            elif xlr > np.real(ck1):
                sign = 1.0
            else:
                dgam = cgam2 - cgam1
                sign = 0.0
                
        if sign != 0.0:
            dgam = 1.0 / (xlam * xlam)
            dgam = sign * ((ct3 * dgam + ct2) * dgam + ct1) / xlam
    
    # Calculate denominators
    den2 = cksm * dgam / (cgam2 * (ck1sq * cgam2 + ck2sq * cgam1))
    den1 = 1.0 / (cgam1 + cgam2) - cksm / cgam2
    
    # Calculate exponential factor
    com = np.exp(-cgam2 * zph)
    
    # Calculate results
    ans = [complex(0.0, 0.0)] * 6
    ans[5] = com * b0 * den1 / ck1
    com = com * den2
    
    if rho == 0.0:
        ans[0] = -com * xlam * xlam * 0.5
        ans[3] = ans[0]
    else:
        b0p = b0p / rho
        ans[0] = -com * xlam * (b0p + b0 * xlam)
        ans[3] = com * xlam * b0p
        
    ans[1] = com * cgam2 * cgam2 * b0
    ans[2] = -ans[3] * cgam2 * rho
    ans[4] = com * b0
    
    return ans

def gshank(start: complex, dela: complex, seed: List[complex], nans: int,
           rho: float, zph: float, ck1: complex, ck2: float, 
           ck1sq: complex, ck2sq: float) -> List[complex]:
    """
    Shanks algorithm for accelerating convergence (GSHANK from Fortran)
    
    Args:
        start: Starting point
        dela: Step size
        seed: Initial values
        nans: Number of answers
        rho, zph: Coordinates
        ck1, ck2: Wave numbers
        ck1sq, ck2sq: Wave numbers squared
        
    Returns:
        Accelerated sum values
    """
    # Initialize arrays
    q1 = [[complex(0.0, 0.0)] * 20 for _ in range(6)]
    q2 = [[complex(0.0, 0.0)] * 20 for _ in range(6)]
    ans1 = [complex(0.0, 0.0)] * 6
    ans2 = seed.copy()
    
    # Integration parameters
    crit = 1e-4
    maxh = 20
    
    # Main integration loop
    b = start
    for int_step in range(maxh):
        a = b
        b = b + dela
        
        # Calculate first step
        sum_vals = rom1_bessel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
        for i in range(nans):
            ans1[i] = ans2[i] + sum_vals[i]
            
        # Calculate second step
        a = b
        b = b + dela
        sum_vals = rom1_bessel(a, b, rho, zph, ck1, ck2, ck1sq, ck2sq)
        for i in range(nans):
            ans2[i] = ans1[i] + sum_vals[i]
            
        # Apply Shanks acceleration
        den = 0.0
        for i in range(nans):
            as1 = ans1[i]
            as2 = ans2[i]
            
            if int_step >= 1:
                for j in range(1, int_step + 1):
                    jm = j - 1
                    aa = q2[i][jm]
                    a1 = q1[i][jm] + as1 - 2.0 * aa
                    if np.real(a1) == 0.0 and np.imag(a1) == 0.0:
                        a1 = q1[i][jm]
                    else:
                        a2 = aa - q1[i][jm]
                        a1 = q1[i][jm] - a2 * a2 / a1
                        
                    a2 = aa + as2 - 2.0 * as1
                    if np.real(a2) == 0.0 and np.imag(a2) == 0.0:
                        a2 = aa
                    else:
                        a2 = aa - (as1 - aa) * (as1 - aa) / a2
                        
                    q1[i][jm] = as1
                    q2[i][jm] = as2
                    as1 = a1
                    as2 = a2
                    
            q1[i][int_step] = as1
            q2[i][int_step] = as2
            amg = abs(np.real(as2)) + abs(np.imag(as2))
            if amg > den:
                den = amg
                
        # Check convergence
        denm = 1e-3 * den * crit
        jm = max(0, int_step - 3)
        converged = True
        
        for j in range(jm, int_step + 1):
            for i in range(nans):
                a1 = q2[i][j]
                den_val = (abs(np.real(a1)) + abs(np.imag(a1))) * crit
                if den_val < denm:
                    den_val = denm
                a1 = q1[i][j] - a1
                amg = abs(np.real(a1)) + abs(np.imag(a1))
                if amg > den_val:
                    converged = False
                    break
            if not converged:
                break
                
        if converged:
            # Return final result
            result = [complex(0.0, 0.0)] * 6
            for i in range(nans):
                result[i] = 0.5 * (q1[i][int_step] + q2[i][int_step])
            return result
    
    # If not converged, return best estimate
    result = [complex(0.0, 0.0)] * 6
    for i in range(nans):
        result[i] = 0.5 * (q1[i][maxh-1] + q2[i][maxh-1])
    return result

# ============================================================================
# BESSEL AND HANKEL FUNCTIONS
# ============================================================================

def bessel(z: complex) -> Tuple[complex, complex]:
    """
    Calculate Bessel function J0 and its derivative (BESSEL from Fortran)
    
    Args:
        z: Complex argument
        
    Returns:
        J0, J0': Bessel function and its derivative
    """
    zms = z * np.conj(z)
    if zms < 1e-12:
        return complex(1.0, 0.0), -0.5 * z
        
    # Series expansion for small arguments
    if zms <= 37.21:
        # Initialize series
        j0 = complex(1.0, 0.0)
        j0p = j0
        zk = j0
        zi = z * z
        
        # Calculate series terms
        for k in range(1, 25):
            zk = zk * (-0.25 / (k * k)) * zi
            j0 = j0 + zk
            j0p = j0p + (1.0 / (k + 1.0)) * zk
            
        j0p = -0.5 * z * j0p
        return j0, j0p
    else:
        # Asymptotic expansion for large arguments
        zi = 1.0 / z
        zi2 = zi * zi
        
        # Polynomial coefficients
        p0z = 1.0 + (0.1121520996 * zi2 - 0.0703125) * zi2
        p1z = 1.0 + (0.1171875 - 0.1441955566 * zi2) * zi2
        q0z = (0.0732421875 * zi2 - 0.125) * zi
        q1z = (0.375 - 0.1025390625 * zi2) * zi
        
        zk = np.exp(1j * (z - 0.7853981635))
        zk = 0.7978845608 * np.sqrt(zi) * zk
        
        j0 = zk * (p0z * np.cos(z - 0.7853981635) - q0z * np.sin(z - 0.7853981635))
        j0p = -zk * (p1z * np.sin(z - 0.7853981635) + q1z * np.cos(z - 0.7853981635))
        
        return j0, j0p

def hankel(z: complex) -> Tuple[complex, complex]:
    """
    Calculate Hankel function H0 and its derivative (HANKEL from Fortran)
    
    Args:
        z: Complex argument
        
    Returns:
        H0, H0': Hankel function and its derivative
    """
    if z == 0.0:
        raise ValueError("Hankel function not valid for z=0")
        
    zms = z * np.conj(z)
    
    # Series expansion for small arguments
    if zms <= 16.81:
        # Initialize series
        j0 = complex(1.0, 0.0)
        j0p = j0
        y0 = complex(0.0, 0.0)
        y0p = y0
        zk = j0
        zi = z * z
        
        # Calculate series terms
        psi = -0.5772156649  # Euler's constant
        for k in range(1, 25):
            zk = zk * (-0.25 / (k * k)) * zi
            j0 = j0 + zk
            j0p = j0p + (1.0 / (k + 1.0)) * zk
            psi = psi + 1.0 / k
            y0 = y0 + (2.0 * psi) * zk
            y0p = y0p + (2.0 * psi + 1.0 / (k + 1.0)) / (k + 1.0) * zk
            
        j0p = -0.5 * z * j0p
        clogz = np.log(0.5 * z)
        y0 = (2.0 * j0 * clogz - y0) / PI + 0.3674669052
        y0p = (2.0 / z + 2.0 * j0p * clogz + 0.5 * y0p * z) / PI - 0.0245785095 * z
        
        h0 = j0 + 1j * y0
        h0p = j0p + 1j * y0p
        return h0, h0p
    else:
        # Asymptotic expansion for large arguments
        zi = 1.0 / z
        zi2 = zi * zi
        
        # Polynomial coefficients
        p0z = 1.0 + (0.1121520996 * zi2 - 0.0703125) * zi2
        p1z = 1.0 + (0.1171875 - 0.1441955566 * zi2) * zi2
        q0z = (0.0732421875 * zi2 - 0.125) * zi
        q1z = (0.375 - 0.1025390625 * zi2) * zi
        
        zk = np.exp(1j * (z - 0.7853981635)) * np.sqrt(zi) * 0.7978845608
        h0 = zk * (p0z + 1j * q0z)
        h0p = 1j * zk * (p1z + 1j * q1z)
        
        return h0, h0p

# ============================================================================
# GROUND REFLECTION COEFFICIENTS
# ============================================================================

def calculate_ground_reflection_coefficients(theta: float, phi: float, 
                                           nec_data: NECData) -> Tuple[complex, complex]:
    """
    Calculate ground reflection coefficients for given angles
    
    Args:
        theta: Elevation angle (radians)
        phi: Azimuth angle (radians)
        nec_data: Main NEC data structure
        
    Returns:
        rrv, rrh: Vertical and horizontal reflection coefficients
    """
    if nec_data.ground.iperf == 1:
        # Perfect ground
        rrv = complex(-1.0, 0.0)
        rrh = complex(-1.0, 0.0)
    else:
        # Finite ground
        cth = np.cos(theta)
        sth = np.sin(theta)
        
        # Calculate complex dielectric constant
        epsc = complex(nec_data.ground.epsr, -nec_data.ground.sig * nec_data.geometry.wlam * 59.96)
        zrati = 1.0 / np.sqrt(epsc)
        
        # Calculate reflection coefficients
        zrsin = np.sqrt(1.0 - zrati * zrati * sth * sth)
        rrv = -(cth - zrati * zrsin) / (cth + zrati * zrsin)
        rrh = (zrati * cth - zrsin) / (zrati * cth + zrsin)
    
    return rrv, rrh

def apply_ground_reflection(ex: complex, ey: complex, ez: complex,
                           rrv: complex, rrh: complex, 
                           theta: float, phi: float) -> Tuple[complex, complex, complex]:
    """
    Apply ground reflection to field components
    
    Args:
        ex, ey, ez: Original field components
        rrv, rrh: Reflection coefficients
        theta, phi: Observation angles
        
    Returns:
        ex_ref, ey_ref, ez_ref: Reflected field components
    """
    # Calculate direction cosines
    cth = np.cos(theta)
    sth = np.sin(theta)
    cph = np.cos(phi)
    sph = np.sin(phi)
    
    # Transform to spherical coordinates
    er = ex * sth * cph + ey * sth * sph + ez * cth
    et = ex * cth * cph + ey * cth * sph - ez * sth
    ep = -ex * sph + ey * cph
    
    # Apply reflection coefficients
    er_ref = rrv * er
    et_ref = rrh * et
    ep_ref = rrh * ep
    
    # Transform back to Cartesian coordinates
    ex_ref = er_ref * sth * cph + et_ref * cth * cph - ep_ref * sph
    ey_ref = er_ref * sth * sph + et_ref * cth * sph + ep_ref * cph
    ez_ref = er_ref * cth - et_ref * sth
    
    return ex_ref, ey_ref, ez_ref

# ============================================================================
# GROUND SCREEN CALCULATIONS
# ============================================================================

def calculate_ground_screen_impedance(rho: float, nec_data: NECData) -> complex:
    """
    Calculate ground screen impedance for radial wire ground screen
    
    Args:
        rho: Radial distance
        nec_data: Main NEC data structure
        
    Returns:
        Ground screen impedance
    """
    if nec_data.ground.nradl == 0:
        return complex(0.0, 0.0)
        
    # Ground screen parameters
    scrwl = nec_data.ground.scrwl
    scrwr = nec_data.ground.scrwr
    nradl = nec_data.ground.nradl
    
    if rho > scrwl:
        return complex(0.0, 0.0)
        
    # Calculate ground screen impedance
    t1 = 1j * 2367.067 / nradl
    t2 = scrwr * nradl
    
    zscrn = t1 * rho * np.log(rho / t2)
    
    # Combine with ground impedance
    zground = complex(nec_data.ground.epsr, -nec_data.ground.sig * nec_data.geometry.wlam * 59.96)
    ztotal = (zscrn * zground) / (ETA * zground + zscrn)
    
    return ztotal

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def test_convergence(f1r: float, f2r: float, f1i: float, f2i: float, dmin: float = 0.0) -> Tuple[float, float]:
    """
    Test for convergence in numerical integration (TEST from Fortran)
    
    Args:
        f1r, f2r: Real parts of function values
        f1i, f2i: Imaginary parts of function values
        dmin: Minimum denominator
        
    Returns:
        tr, ti: Relative errors
    """
    den = abs(f2r)
    tr = abs(f2i)
    if den < tr:
        den = tr
    if den < dmin:
        den = dmin
    if den < 1e-37:
        return 0.0, 0.0
        
    tr = abs((f1r - f2r) / den)
    ti = abs((f1i - f2i) / den)
    
    return tr, ti

# Global variables for Sommerfeld calculations
ct1 = complex(0.0, 0.0)
ct2 = complex(0.0, 0.0)
ct3 = complex(0.0, 0.0)
cksm = complex(0.0, 0.0) 