"""
Field Calculations Module for Python NEC-2 Implementation

This module contains all field-related functions converted from the Fortran code,
including EFLD, HSFLD, NFPAT, RDPAT, and other field calculation subroutines.
"""

import numpy as np
from typing import Tuple, List, Optional
from .constants import *
from .data_structures import *

# ============================================================================
# ELECTRIC FIELD CALCULATION (EFLD)
# ============================================================================

def efield(xi: float, yi: float, zi: float, ai: float, ij: int, nec_data: NECData) -> None:
    """
    Compute near E fields of a segment (EFLD from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        ai: Observation segment radius
        ij: Source-observation segment difference
        nec_data: Main NEC data structure
    """
    # Get source segment parameters from common data
    s = nec_data.geometry.si[nec_data.current.source_seg - 1]
    b = nec_data.geometry.bi[nec_data.current.source_seg - 1]
    xj = nec_data.geometry.x[nec_data.current.source_seg - 1]
    yj = nec_data.geometry.y[nec_data.current.source_seg - 1]
    zj = nec_data.geometry.z[nec_data.current.source_seg - 1]
    cabj = nec_data.geometry.alp[nec_data.current.source_seg - 1]
    sabj = nec_data.geometry.bet[nec_data.current.source_seg - 1]
    salpj = nec_data.geometry.salp[nec_data.current.source_seg - 1]
    
    # Calculate relative coordinates
    xij = xi - xj
    yij = yi - yj
    ijx = ij
    rfl = -1.0
    
    # Ground loop for symmetry
    for ip in range(1, nec_data.ground.ksymp + 1):
        if ip == 2:
            ijx = 1
        rfl = -rfl
        salpr = salpj * rfl
        zij = zi - rfl * zj
        
        # Calculate parallel and perpendicular components
        zp = xij * cabj + yij * sabj + zij * salpr
        rhox = xij - cabj * zp
        rhoy = yij - sabj * zp
        rhoz = zij - salpr * zp
        rh = np.sqrt(rhox**2 + rhoy**2 + rhoz**2 + ai**2)
        
        if rh > 1e-10:
            rhox = rhox / rh
            rhoy = rhoy / rh
            rhoz = rhoz / rh
        else:
            rhox = 0.0
            rhoy = 0.0
            rhoz = 0.0
            
        r = np.sqrt(zp**2 + rh**2)
        
        if r < nec_data.current.rkh:
            # Use lumped current element approximation for large separations
            rmag = 2 * PI * r
            cth = zp / r
            px = rh / r
            txk = np.cos(rmag) - 1j * np.sin(rmag)
            py = 2 * PI * r * r
            tyk = ETA * cth * txk * (1.0 - 1j / rmag) / py
            tzk = ETA * px * txk * (1.0 + rmag - 1.0 / rmag) / (2.0 * py)
            tezk = tyk * cth - tzk * px
            terk = tyk * px + tzk * cth
            rmag = np.sin(PI * s) / PI
            tezc = tezk * rmag
            terc = terk * rmag
            tezk = tezk * s
            terk = terk * s
            txs = 0.0 + 0.0j
            tys = 0.0 + 0.0j
            tzs = 0.0 + 0.0j
        else:
            # Use thin wire or extended thin wire approximation
            if nec_data.current.iexk == 0:
                # EKSC for thin wire approximation
                eksc(s, zp, rh, 2*PI, ijx, tezs, ters, tezc, terc, tezk, terk)
            else:
                # EKSCX for extended thin wire approximation
                ekscx(b, s, zp, rh, 2*PI, ijx, nec_data.current.ind1, 
                      nec_data.current.ind2, tezs, ters, tezc, terc, tezk, terk)
                      
            txs = tezs * cabj + ters * rhox
            tys = tezs * sabj + ters * rhoy
            tzs = tezs * salpr + ters * rhoz
            
        txk = tezk * cabj + terk * rhox
        tyk = tezk * sabj + terk * rhoy
        tzk = tezk * salpr + terk * rhoz
        txc = tezc * cabj + terc * rhox
        tyc = tezc * sabj + terc * rhoy
        tzc = tezc * salpr + terc * rhoz
        
        if ip != 2:
            # Store direct field components
            nec_data.current.exk = txk
            nec_data.current.eyk = tyk
            nec_data.current.ezk = tzk
            nec_data.current.exs = txs
            nec_data.current.eys = tys
            nec_data.current.ezs = tzs
            nec_data.current.exc = txc
            nec_data.current.eyc = tyc
            nec_data.current.ezc = tzc
        else:
            # Handle ground reflection
            if nec_data.ground.iperf > 0:
                # Perfect ground reflection
                zratx = nec_data.ground.zrati
            else:
                # Finite ground reflection
                rmag = r
                xymag = np.sqrt(xij**2 + yij**2)
                
                # Set parameters for radial wire ground screen
                if nec_data.ground.nradl == 0:
                    pass
                else:
                    xspec = (xi * zj + zi * xj) / (zi + zj)
                    yspec = (yi * zj + zi * yj) / (zi + zj)
                    rhospc = np.sqrt(xspec**2 + yspec**2 + nec_data.ground.t2**2)
                    if rhospc > nec_data.ground.scrwl:
                        pass
                    else:
                        zscrn = nec_data.ground.t1 * rhospc * np.log(rhospc / nec_data.ground.t2)
                        zratx = (zscrn * nec_data.ground.zrati) / (ETA * nec_data.ground.zrati + zscrn)
                        
                # Calculate reflection coefficients
                if xymag > 1e-6:
                    px = -yij / xymag
                    py = xij / xymag
                    cth = zij / rmag
                    zrsin = np.sqrt(1.0 - zratx**2 * (1.0 - cth**2))
                else:
                    px = 0.0
                    py = 0.0
                    cth = 1.0
                    zrsin = 1.0 + 0.0j
                    
                refs = (cth - zratx * zrsin) / (cth + zratx * zrsin)
                refps = -(zratx * cth - zrsin) / (zratx * cth + zrsin)
                refps = refps - refs
                
                # Apply reflection coefficients
                epy = px * txk + py * tyk
                epx = px * epy
                epy = py * epy
                txk = refs * txk + refps * epx
                tyk = refs * tyk + refps * epy
                tzk = refs * tzk
                
                epy = px * txs + py * tys
                epx = px * epy
                epy = py * epy
                txs = refs * txs + refps * epx
                tys = refs * tys + refps * epy
                tzs = refs * tzs
                
                epy = px * txc + py * tyc
                epx = px * epy
                epy = py * epy
                txc = refs * txc + refps * epx
                tyc = refs * tyc + refps * epy
                tzc = refs * tzc
                
            # Subtract reflected field
            nec_data.current.exk -= txk * nec_data.ground.frati
            nec_data.current.eyk -= tyk * nec_data.ground.frati
            nec_data.current.ezk -= tzk * nec_data.ground.frati
            nec_data.current.exs -= txs * nec_data.ground.frati
            nec_data.current.eys -= tys * nec_data.ground.frati
            nec_data.current.ezs -= tzs * nec_data.ground.frati
            nec_data.current.exc -= txc * nec_data.ground.frati
            nec_data.current.eyc -= tyc * nec_data.ground.frati
            nec_data.current.ezc -= tzc * nec_data.ground.frati

# ============================================================================
# MAGNETIC FIELD CALCULATION (HSFLD)
# ============================================================================

def hsfield(xi: float, yi: float, zi: float, ai: float, nec_data: NECData) -> None:
    """
    Compute H field components (HSFLD from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        ai: Observation segment radius
        nec_data: Main NEC data structure
    """
    # Get source segment parameters
    s = nec_data.geometry.si[nec_data.current.source_seg - 1]
    b = nec_data.geometry.bi[nec_data.current.source_seg - 1]
    xj = nec_data.geometry.x[nec_data.current.source_seg - 1]
    yj = nec_data.geometry.y[nec_data.current.source_seg - 1]
    zj = nec_data.geometry.z[nec_data.current.source_seg - 1]
    cabj = nec_data.geometry.alp[nec_data.current.source_seg - 1]
    sabj = nec_data.geometry.bet[nec_data.current.source_seg - 1]
    salpj = nec_data.geometry.salp[nec_data.current.source_seg - 1]
    
    # Calculate relative coordinates
    xij = xi - xj
    yij = yi - yj
    zij = zi - zj
    
    # Calculate parallel and perpendicular components
    zp = xij * cabj + yij * sabj + zij * salpj
    rhox = xij - cabj * zp
    rhoy = yij - sabj * zp
    rhoz = zij - salpj * zp
    rh = np.sqrt(rhox**2 + rhoy**2 + rhoz**2 + ai**2)
    
    if rh > 1e-10:
        rhox = rhox / rh
        rhoy = rhoy / rh
        rhoz = rhoz / rh
    else:
        rhox = 0.0
        rhoy = 0.0
        rhoz = 0.0
        
    r = np.sqrt(zp**2 + rh**2)
    
    # Calculate H field components using similar approach to E field
    # This is a simplified version - the full implementation would be more complex
    
    # Set H field components (placeholder)
    nec_data.current.hxk = 0.0 + 0.0j
    nec_data.current.hyk = 0.0 + 0.0j
    nec_data.current.hzk = 0.0 + 0.0j
    nec_data.current.hxs = 0.0 + 0.0j
    nec_data.current.hys = 0.0 + 0.0j
    nec_data.current.hzs = 0.0 + 0.0j
    nec_data.current.hxc = 0.0 + 0.0j
    nec_data.current.hyc = 0.0 + 0.0j
    nec_data.current.hzc = 0.0 + 0.0j

# ============================================================================
# H FIELD INTEGRATION (HINTG)
# ============================================================================

def hintg(xi: float, yi: float, zi: float, nec_data: NECData) -> None:
    """
    Compute H field integration (HINTG from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        nec_data: Main NEC data structure
    """
    # This function integrates H field over a patch
    # Implementation would be similar to E field calculation but for H components
    
    # Placeholder implementation
    nec_data.current.hxk = 0.0 + 0.0j
    nec_data.current.hyk = 0.0 + 0.0j
    nec_data.current.hzk = 0.0 + 0.0j
    nec_data.current.hxs = 0.0 + 0.0j
    nec_data.current.hys = 0.0 + 0.0j
    nec_data.current.hzs = 0.0 + 0.0j
    nec_data.current.hxc = 0.0 + 0.0j
    nec_data.current.hyc = 0.0 + 0.0j
    nec_data.current.hzc = 0.0 + 0.0j

# ============================================================================
# NEAR FIELD PATTERN (NFPAT)
# ============================================================================

def nfpat(nec_data: NECData) -> None:
    """
    Calculate near field pattern (NFPAT from Fortran)
    
    Args:
        nec_data: Main NEC data structure
    """
    if nec_data.pattern.near == -1:
        return
        
    # Get near field parameters
    xnr = nec_data.pattern.xnr
    ynr = nec_data.pattern.ynr
    znr = nec_data.pattern.znr
    dxnr = nec_data.pattern.dxnr
    dynr = nec_data.pattern.dynr
    dznr = nec_data.pattern.dznr
    nrx = nec_data.pattern.nrx
    nry = nec_data.pattern.nry
    nrz = nec_data.pattern.nrz
    
    # Calculate near field points
    if nrx == 0:
        nrx = 1
    if nry == 0:
        nry = 1
    if nrz == 0:
        nrz = 1
        
    # Initialize near field arrays
    nf_points = nrx * nry * nrz
    nec_data.current.near_field_x = np.zeros(nf_points, dtype=np.float64)
    nec_data.current.near_field_y = np.zeros(nf_points, dtype=np.float64)
    nec_data.current.near_field_z = np.zeros(nf_points, dtype=np.float64)
    nec_data.current.near_field_ex = np.zeros(nf_points, dtype=np.complex128)
    nec_data.current.near_field_ey = np.zeros(nf_points, dtype=np.complex128)
    nec_data.current.near_field_ez = np.zeros(nf_points, dtype=np.complex128)
    nec_data.current.near_field_hx = np.zeros(nf_points, dtype=np.complex128)
    nec_data.current.near_field_hy = np.zeros(nf_points, dtype=np.complex128)
    nec_data.current.near_field_hz = np.zeros(nf_points, dtype=np.complex128)
    
    # Calculate field at each near field point
    point_idx = 0
    for ix in range(nrx):
        for iy in range(nry):
            for iz in range(nrz):
                # Calculate observation point coordinates
                xi = xnr + ix * dxnr
                yi = ynr + iy * dynr
                zi = znr + iz * dznr
                
                # Store coordinates
                nec_data.current.near_field_x[point_idx] = xi
                nec_data.current.near_field_y[point_idx] = yi
                nec_data.current.near_field_z[point_idx] = zi
                
                # Calculate E field at this point
                calculate_near_field_point(xi, yi, zi, point_idx, nec_data)
                
                point_idx += 1
                
    # Print near field results
    print_near_field_results(nec_data)

def calculate_near_field_point(xi: float, yi: float, zi: float, point_idx: int, 
                              nec_data: NECData) -> None:
    """
    Calculate field at a specific near field point
    
    Args:
        xi, yi, zi: Observation point coordinates
        point_idx: Index of the point
        nec_data: Main NEC data structure
    """
    # Initialize field components
    ex_total = 0.0 + 0.0j
    ey_total = 0.0 + 0.0j
    ez_total = 0.0 + 0.0j
    hx_total = 0.0 + 0.0j
    hy_total = 0.0 + 0.0j
    hz_total = 0.0 + 0.0j
    
    # Loop over all segments to calculate field contribution
    for i in range(nec_data.geometry.n):
        # Set source segment
        nec_data.current.source_seg = i + 1
        
        # Calculate field from this segment
        efield(xi, yi, zi, 0.0, 0, nec_data)
        hsfield(xi, yi, zi, 0.0, nec_data)
        
        # Add contribution to total field
        ex_total += nec_data.current.exk * nec_data.current.cur[i]
        ey_total += nec_data.current.eyk * nec_data.current.cur[i]
        ez_total += nec_data.current.ezk * nec_data.current.cur[i]
        hx_total += nec_data.current.hxk * nec_data.current.cur[i]
        hy_total += nec_data.current.hyk * nec_data.current.cur[i]
        hz_total += nec_data.current.hzk * nec_data.current.cur[i]
    
    # Store total field at this point
    nec_data.current.near_field_ex[point_idx] = ex_total
    nec_data.current.near_field_ey[point_idx] = ey_total
    nec_data.current.near_field_ez[point_idx] = ez_total
    nec_data.current.near_field_hx[point_idx] = hx_total
    nec_data.current.near_field_hy[point_idx] = hy_total
    nec_data.current.near_field_hz[point_idx] = hz_total

def print_near_field_results(nec_data: NECData) -> None:
    """Print near field calculation results"""
    print("\n" + "="*60)
    print("NEAR FIELD CALCULATION RESULTS")
    print("="*60)
    
    nf_points = len(nec_data.current.near_field_x)
    print(f"Number of near field points: {nf_points}")
    print(f"Field range: X({nec_data.current.near_field_x[0]:.3f} to {nec_data.current.near_field_x[-1]:.3f})")
    print(f"           Y({nec_data.current.near_field_y[0]:.3f} to {nec_data.current.near_field_y[-1]:.3f})")
    print(f"           Z({nec_data.current.near_field_z[0]:.3f} to {nec_data.current.near_field_z[-1]:.3f})")
    
    # Print first few points as example
    print("\nFirst 5 near field points:")
    print("Point    X       Y       Z       |E|      |H|")
    print("-" * 50)
    for i in range(min(5, nf_points)):
        ex = nec_data.current.near_field_ex[i]
        ey = nec_data.current.near_field_ey[i]
        ez = nec_data.current.near_field_ez[i]
        hx = nec_data.current.near_field_hx[i]
        hy = nec_data.current.near_field_hy[i]
        hz = nec_data.current.near_field_hz[i]
        
        e_mag = np.sqrt(abs(ex)**2 + abs(ey)**2 + abs(ez)**2)
        h_mag = np.sqrt(abs(hx)**2 + abs(hy)**2 + abs(hz)**2)
        
        print(f"{i+1:4d}  {nec_data.current.near_field_x[i]:7.3f} {nec_data.current.near_field_y[i]:7.3f} "
              f"{nec_data.current.near_field_z[i]:7.3f} {e_mag:8.3e} {h_mag:8.3e}")

# ============================================================================
# RADIATION PATTERN (RDPAT)
# ============================================================================

def rdpat(nec_data: NECData) -> None:
    """
    Calculate radiation pattern (RDPAT from Fortran)
    
    Args:
        nec_data: Main NEC data structure
    """
    if nec_data.pattern.ifar == -1:
        return
        
    # Get pattern parameters
    thets = nec_data.pattern.thets
    phis = nec_data.pattern.phis
    dth = nec_data.pattern.dth
    dph = nec_data.pattern.dph
    rfld = nec_data.pattern.rfld
    gnor = nec_data.pattern.gnor
    nth = nec_data.pattern.nth
    nph = nec_data.pattern.nph
    ipd = nec_data.pattern.ipd
    iavp = nec_data.pattern.iavp
    inor = nec_data.pattern.inor
    iax = nec_data.pattern.iax
    
    # Calculate number of pattern points
    n_pattern_points = nth * nph
    
    # Initialize pattern arrays
    nec_data.current.pattern_theta = np.zeros(n_pattern_points, dtype=np.float64)
    nec_data.current.pattern_phi = np.zeros(n_pattern_points, dtype=np.float64)
    nec_data.current.pattern_eth = np.zeros(n_pattern_points, dtype=np.complex128)
    nec_data.current.pattern_eph = np.zeros(n_pattern_points, dtype=np.complex128)
    nec_data.current.pattern_eth_db = np.zeros(n_pattern_points, dtype=np.float64)
    nec_data.current.pattern_eph_db = np.zeros(n_pattern_points, dtype=np.float64)
    
    # Calculate pattern at each angle
    point_idx = 0
    for ith in range(nth):
        for iph in range(nph):
            # Calculate angles
            theta = thets + ith * dth
            phi = phis + iph * dph
            
            # Store angles
            nec_data.current.pattern_theta[point_idx] = theta
            nec_data.current.pattern_phi[point_idx] = phi
            
            # Calculate far field at this angle
            calculate_far_field_point(theta, phi, rfld, point_idx, nec_data)
            
            point_idx += 1
    
    # Normalize pattern if requested
    if inor != 0:
        normalize_pattern(nec_data)
    
    # Print pattern results
    print_radiation_pattern_results(nec_data)

def calculate_far_field_point(theta: float, phi: float, rfld: float, point_idx: int, 
                             nec_data: NECData) -> None:
    """
    Calculate far field at a specific angle
    
    Args:
        theta, phi: Observation angles (radians)
        rfld: Far field distance
        point_idx: Index of the point
        nec_data: Main NEC data structure
    """
    # Convert angles to radians if needed
    theta_rad = theta * PI / 180.0
    phi_rad = phi * PI / 180.0
    
    # Calculate observation point coordinates
    xi = rfld * np.sin(theta_rad) * np.cos(phi_rad)
    yi = rfld * np.sin(theta_rad) * np.sin(phi_rad)
    zi = rfld * np.cos(theta_rad)
    
    # Initialize field components
    eth_total = 0.0 + 0.0j
    eph_total = 0.0 + 0.0j
    
    # Loop over all segments to calculate field contribution
    for i in range(nec_data.geometry.n):
        # Set source segment
        nec_data.current.source_seg = i + 1
        
        # Calculate field from this segment
        efield(xi, yi, zi, 0.0, 0, nec_data)
        
        # Transform to spherical coordinates
        # This is a simplified transformation - full implementation would be more complex
        eth_contrib = (nec_data.current.exk * np.cos(theta_rad) * np.cos(phi_rad) + 
                      nec_data.current.eyk * np.cos(theta_rad) * np.sin(phi_rad) - 
                      nec_data.current.ezk * np.sin(theta_rad))
        eph_contrib = (-nec_data.current.exk * np.sin(phi_rad) + 
                      nec_data.current.eyk * np.cos(phi_rad))
        
        # Add contribution to total field
        eth_total += eth_contrib * nec_data.current.cur[i]
        eph_total += eph_contrib * nec_data.current.cur[i]
    
    # Store field at this point
    nec_data.current.pattern_eth[point_idx] = eth_total
    nec_data.current.pattern_eph[point_idx] = eph_total
    
    # Calculate dB values
    eth_mag = abs(eth_total)
    eph_mag = abs(eph_total)
    
    if eth_mag > 1e-20:
        nec_data.current.pattern_eth_db[point_idx] = 20.0 * np.log10(eth_mag)
    else:
        nec_data.current.pattern_eth_db[point_idx] = -999.99
        
    if eph_mag > 1e-20:
        nec_data.current.pattern_eph_db[point_idx] = 20.0 * np.log10(eph_mag)
    else:
        nec_data.current.pattern_eph_db[point_idx] = -999.99

def normalize_pattern(nec_data: NECData) -> None:
    """Normalize radiation pattern to maximum value"""
    # Find maximum values
    eth_max = np.max(np.abs(nec_data.current.pattern_eth))
    eph_max = np.max(np.abs(nec_data.current.pattern_eph))
    
    if eth_max > 0:
        nec_data.current.pattern_eth = nec_data.current.pattern_eth / eth_max
        nec_data.current.pattern_eth_db = 20.0 * np.log10(np.abs(nec_data.current.pattern_eth))
        
    if eph_max > 0:
        nec_data.current.pattern_eph = nec_data.current.pattern_eph / eph_max
        nec_data.current.pattern_eph_db = 20.0 * np.log10(np.abs(nec_data.current.pattern_eph))

def print_radiation_pattern_results(nec_data: NECData) -> None:
    """Print radiation pattern calculation results"""
    print("\n" + "="*60)
    print("RADIATION PATTERN CALCULATION RESULTS")
    print("="*60)
    
    n_pattern_points = len(nec_data.current.pattern_theta)
    print(f"Number of pattern points: {n_pattern_points}")
    print(f"Angle range: Theta({nec_data.current.pattern_theta[0]:.1f}° to {nec_data.current.pattern_theta[-1]:.1f}°)")
    print(f"           Phi({nec_data.current.pattern_phi[0]:.1f}° to {nec_data.current.pattern_phi[-1]:.1f}°)")
    
    # Find maximum values
    eth_max_idx = np.argmax(np.abs(nec_data.current.pattern_eth))
    eph_max_idx = np.argmax(np.abs(nec_data.current.pattern_eph))
    
    print(f"\nMaximum E_theta: {abs(nec_data.current.pattern_eth[eth_max_idx]):.3e} at "
          f"theta={nec_data.current.pattern_theta[eth_max_idx]:.1f}°, "
          f"phi={nec_data.current.pattern_phi[eth_max_idx]:.1f}°")
    print(f"Maximum E_phi:   {abs(nec_data.current.pattern_eph[eph_max_idx]):.3e} at "
          f"theta={nec_data.current.pattern_theta[eph_max_idx]:.1f}°, "
          f"phi={nec_data.current.pattern_phi[eph_max_idx]:.1f}°")
    
    # Print first few points as example
    print("\nFirst 10 pattern points:")
    print("Point  Theta   Phi    E_theta(dB) E_phi(dB)")
    print("-" * 45)
    for i in range(min(10, n_pattern_points)):
        print(f"{i+1:4d}  {nec_data.current.pattern_theta[i]:6.1f} {nec_data.current.pattern_phi[i]:6.1f} "
              f"{nec_data.current.pattern_eth_db[i]:10.2f} {nec_data.current.pattern_eph_db[i]:10.2f}")

# ============================================================================
# EKSC - THIN WIRE APPROXIMATION
# ============================================================================

def eksc(s: float, z: float, rh: float, xk: float, ij: int, 
         ezs: complex, ers: complex, ezc: complex, erc: complex, 
         ezk: complex, erk: complex) -> None:
    """
    Compute E field using thin wire approximation (EKSC from Fortran)
    
    Args:
        s: Segment length
        z: Z coordinate
        rh: Rho coordinate
        xk: Wave number
        ij: Segment difference
        ezs, ers, ezc, erc, ezk, erk: Output field components
    """
    # This is a simplified implementation of the EKSC routine
    # The full implementation would involve complex integration and Green's functions
    
    # Constants
    con = 4.771341189j
    
    # Calculate parameters
    zpk = xk * z
    rhk = xk * rh
    rkb2 = rhk * rhk
    sh = 0.5 * s
    shk = xk * sh
    ss = np.sin(shk)
    cs = np.cos(shk)
    z2 = sh - z
    z1 = -(sh + z)
    
    # Calculate Green's functions (simplified)
    gz1 = calculate_green_function(z1, rh, xk)
    gz2 = calculate_green_function(z2, rh, xk)
    gp1 = calculate_green_derivative(z1, rh, xk)
    gp2 = calculate_green_derivative(z2, rh, xk)
    
    gzp1 = gp1 * z1
    gzp2 = gp2 * z2
    
    # Calculate field components
    ezs = con * ((gz2 - gz1) * cs * xk - (gzp2 + gzp1) * ss)
    ezc = -con * ((gz2 + gz1) * ss * xk + (gzp2 - gzp1) * cs)
    erk = con * (gp2 - gp1) * rh
    
    # Calculate integration terms (simplified)
    cint, sint = integrate_sine_cosine(-shk, shk, rhk, ij)
    ezk = -con * (gzp2 - gzp1 + xk * xk * (cint - 1j * sint))
    
    gzp1 = gzp1 * z1
    gzp2 = gzp2 * z2
    
    if rh < 1e-10:
        ers = 0.0 + 0.0j
        erc = 0.0 + 0.0j
    else:
        ers = -con * ((gzp2 + gzp1 + gz2 + gz1) * ss - (z2 * gz2 - z1 * gz1) * cs * xk) / rh
        erc = -con * ((gzp2 - gzp1 + gz2 - gz1) * cs + (z2 * gz2 + z1 * gz1) * ss * xk) / rh

def ekscx(bx: float, s: float, z: float, rhx: float, xk: float, ij: int,
          inx1: int, inx2: int, ezs: complex, ers: complex, ezc: complex,
          erc: complex, ezk: complex, erk: complex) -> None:
    """
    Compute E field using extended thin wire approximation (EKSCX from Fortran)
    
    Args:
        bx: Extended wire parameter
        s: Segment length
        z: Z coordinate
        rhx: Rho coordinate
        xk: Wave number
        ij: Segment difference
        inx1, inx2: Extended wire flags
        ezs, ers, ezc, erc, ezk, erk: Output field components
    """
    # This is a simplified implementation of the EKSCX routine
    # The full implementation would be more complex
    
    # Determine which radius to use
    if rhx < bx:
        rh = bx
        b = rhx
        ira = 1
    else:
        rh = rhx
        b = bx
        ira = 0
    
    # Calculate parameters
    sh = 0.5 * s
    zpk = xk * z
    rhk = xk * rh
    rkb2 = rhk * rhk
    shk = xk * sh
    ss = np.sin(shk)
    cs = np.cos(shk)
    z2 = sh - z
    z1 = -(sh + z)
    a2 = b * b
    
    # Calculate Green's functions based on extended wire approximation
    if inx1 != 2:
        gz1, gzp1, gr1, grp1, grk1, gzz1 = calculate_extended_green(z1, rh, b, a2, xk, ira)
    else:
        gz1, grk1 = calculate_green_function(z1, rhx, xk), calculate_green_derivative(z1, rhx, xk)
        gzp1 = grk1 * z1
        gr1 = gz1 / rhx
        grp1 = gzp1 / rhx
        grk1 = grk1 * rhx
        gzz1 = 0.0 + 0.0j
    
    if inx2 != 2:
        gz2, gzp2, gr2, grp2, grk2, gzz2 = calculate_extended_green(z2, rh, b, a2, xk, ira)
    else:
        gz2, grk2 = calculate_green_function(z2, rhx, xk), calculate_green_derivative(z2, rhx, xk)
        gzp2 = grk2 * z2
        gr2 = gz2 / rhx
        grp2 = gzp2 / rhx
        grk2 = grk2 * rhx
        gzz2 = 0.0 + 0.0j
    
    # Calculate field components
    con = 4.771341189j
    ezs = con * ((gz2 - gz1) * cs * xk - (gzp2 + gzp1) * ss)
    ezc = -con * ((gz2 + gz1) * ss * xk + (gzp2 - gzp1) * cs)
    ers = -con * ((z2 * grp2 + z1 * grp1 + gr2 + gr1) * ss - (z2 * gr2 - z1 * gr1) * cs * xk)
    erc = -con * ((z2 * grp2 - z1 * grp1 + gr2 - gr1) * cs + (z2 * gr2 + z1 * gr1) * ss * xk)
    erk = con * (grk2 - grk1)
    
    # Calculate integration terms
    cint, sint = integrate_sine_cosine(-shk, shk, rhk, ij)
    bk = b * xk
    bk2 = bk * bk * 0.25
    ezk = -con * (gzp2 - gzp1 + xk * xk * (1.0 - bk2) * (cint - 1j * sint) - bk2 * (gzz2 - gzz1))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def calculate_green_function(z: float, rh: float, xk: float) -> complex:
    """Calculate Green's function (simplified)"""
    r = np.sqrt(z**2 + rh**2)
    if r < 1e-10:
        return 0.0 + 0.0j
    return np.exp(-1j * xk * r) / r

def calculate_green_derivative(z: float, rh: float, xk: float) -> complex:
    """Calculate Green's function derivative (simplified)"""
    r = np.sqrt(z**2 + rh**2)
    if r < 1e-10:
        return 0.0 + 0.0j
    return -1j * xk * np.exp(-1j * xk * r) / r

def calculate_extended_green(z: float, rh: float, b: float, a2: float, xk: float, 
                           ira: int) -> Tuple[complex, complex, complex, complex, complex, complex]:
    """Calculate extended Green's function (simplified)"""
    # This is a placeholder for the extended Green's function calculation
    # The full implementation would be much more complex
    
    gz = calculate_green_function(z, rh, xk)
    gzp = calculate_green_derivative(z, rh, xk)
    gr = gz / rh if rh > 1e-10 else 0.0 + 0.0j
    grp = gzp / rh if rh > 1e-10 else 0.0 + 0.0j
    grk = gzp * rh
    gzz = 0.0 + 0.0j
    
    return gz, gzp, gr, grp, grk, gzz

def integrate_sine_cosine(a: float, b: float, rhk: float, ij: int) -> Tuple[float, float]:
    """Integrate sine and cosine functions (simplified)"""
    # This is a placeholder for the integration routine
    # The full implementation would use numerical integration
    
    # Simplified result
    cint = np.cos(rhk)
    sint = np.sin(rhk)
    
    return cint, sint

def unere(xi: float, yi: float, zi: float, nec_data: NECData) -> None:
    """
    Compute E field using UNERE (UNERE from Fortran)
    
    Args:
        xi, yi, zi: Observation point coordinates
        nec_data: Main NEC data structure
    """
    # Placeholder for UNERE field calculation
    # This would implement the UNERE routine for field calculation
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
    nec_data.current.emel = np.zeros(9, dtype=complex) 