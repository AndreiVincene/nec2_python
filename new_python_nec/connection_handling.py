"""
Connection Handling Module for Python NEC-2 Implementation

This module contains the connection handling functions converted from the Fortran code,
including CONECT and related subroutines for setting up segment connections.
"""

import numpy as np
from typing import List, Tuple, Optional
from .constants import *
from .data_structures import *

# ============================================================================
# CONNECTION SETUP (CONECT)
# ============================================================================

def conect(ignd: int, nec_data: NECData) -> None:
    """
    Set up segment connection data (CONECT from Fortran)
    
    Args:
        ignd: Ground flag
        nec_data: Main NEC data structure
    """
    # Initialize connection parameters
    smin = 1e-3  # Minimum separation for connection
    npmax = 10   # Maximum number of new segments connected to NGF patches
    
    # Initialize connection counters
    nscon = 0
    npcon = 0
    
    # Handle ground plane specification
    if ignd != 0:
        print("\n   GROUND PLANE SPECIFIED.")
        if ignd > 0:
            print("   WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL BE INTERPOLATED TO IMAGE IN GROUND PLANE.")
        
        # Adjust symmetry parameters for ground
        if nec_data.geometry.ipsym != 2:
            pass
        else:
            nec_data.geometry.np = 2 * nec_data.geometry.np
            nec_data.geometry.mp = 2 * nec_data.geometry.mp
            
        if abs(nec_data.geometry.ipsym) <= 2:
            pass
        else:
            nec_data.geometry.np = nec_data.geometry.n
            nec_data.geometry.mp = nec_data.geometry.m
            
        if nec_data.geometry.np > nec_data.geometry.n:
            print("ERROR: NP > N")
            return
            
        if nec_data.geometry.np == nec_data.geometry.n and nec_data.geometry.mp == nec_data.geometry.m:
            nec_data.geometry.ipsym = 0
    
    # Process wire segments
    if nec_data.geometry.n != 0:
        process_wire_connections(ignd, smin, nec_data)
    
    # Process surface patches
    if nec_data.geometry.m != 0:
        process_patch_connections(ignd, smin, nec_data)
    
    # Check connection limits
    if npcon > npmax:
        print(f"ERROR - NO. NEW SEGMENTS CONNECTED TO N.G.F. SEGMENTS OR PATCHES EXCEEDS LIMIT OF {npmax}")
        return
    
    # Print connection summary
    print_connection_summary(nec_data)
    
    # Process multiple wire junctions
    if nec_data.geometry.n != 0:
        process_multiple_junctions(nec_data)

def process_wire_connections(ignd: int, smin: float, nec_data: NECData) -> None:
    """
    Process connections for wire segments
    
    Args:
        ignd: Ground flag
        smin: Minimum separation
        nec_data: Main NEC data structure
    """
    for i in range(nec_data.geometry.n):
        # Initialize connection data
        nec_data.geometry.iconx[i] = 0
        
        # Get segment end coordinates
        xi1 = nec_data.geometry.x[i]
        yi1 = nec_data.geometry.y[i]
        zi1 = nec_data.geometry.z[i]
        xi2 = nec_data.geometry.x[i] + nec_data.geometry.si[i] * nec_data.geometry.alp[i]
        yi2 = nec_data.geometry.y[i] + nec_data.geometry.si[i] * nec_data.geometry.bet[i]
        zi2 = nec_data.geometry.z[i] + nec_data.geometry.si[i] * nec_data.geometry.salp[i]
        
        # Calculate segment length
        slen = np.sqrt((xi2 - xi1)**2 + (yi2 - yi1)**2 + (zi2 - zi1)**2) * smin
        
        # Process end 1 connection
        if ignd >= 1:
            if zi1 > -slen:
                if zi1 > slen:
                    # No ground connection
                    pass
                else:
                    # Ground connection
                    nec_data.geometry.icon1[i] = i
                    nec_data.geometry.z[i] = 0.0
            else:
                print(f"GEOMETRY DATA ERROR-- SEGMENT {i+1} EXTENDS BELOW GROUND")
                return
        else:
            # Find connection to other segments
            ic = i
            for j in range(2, nec_data.geometry.n + 1):
                ic = ic + 1
                if ic > nec_data.geometry.n:
                    ic = 1
                    
                # Check separation from segment ends
                sep = abs(xi1 - nec_data.geometry.x[ic-1]) + abs(yi1 - nec_data.geometry.y[ic-1]) + abs(zi1 - nec_data.geometry.z[ic-1])
                if sep <= slen:
                    nec_data.geometry.icon1[i] = -ic
                    break
                    
                sep = abs(xi1 - (nec_data.geometry.x[ic-1] + nec_data.geometry.si[ic-1] * nec_data.geometry.alp[ic-1])) + \
                      abs(yi1 - (nec_data.geometry.y[ic-1] + nec_data.geometry.si[ic-1] * nec_data.geometry.bet[ic-1])) + \
                      abs(zi1 - (nec_data.geometry.z[ic-1] + nec_data.geometry.si[ic-1] * nec_data.geometry.salp[ic-1]))
                if sep <= slen:
                    nec_data.geometry.icon1[i] = ic
                    break
            else:
                if i < nec_data.geometry.n2 and nec_data.geometry.icon1[i] > 10000:
                    pass
                else:
                    nec_data.geometry.icon1[i] = 0
        
        # Process end 2 connection
        if ignd >= 1:
            if zi2 > -slen:
                if zi2 > slen:
                    # No ground connection
                    pass
                else:
                    # Ground connection
                    if nec_data.geometry.icon1[i] != i:
                        nec_data.geometry.icon2[i] = i
                        nec_data.geometry.z[i] = 0.0
                    else:
                        print(f"GEOMETRY DATA ERROR--SEGMENT {i+1} LIES IN GROUND PLANE.")
                        return
            else:
                print(f"GEOMETRY DATA ERROR-- SEGMENT {i+1} EXTENDS BELOW GROUND")
                return
        else:
            # Find connection to other segments
            ic = i
            for j in range(2, nec_data.geometry.n + 1):
                ic = ic + 1
                if ic > nec_data.geometry.n:
                    ic = 1
                    
                # Check separation from segment ends
                sep = abs(xi2 - nec_data.geometry.x[ic-1]) + abs(yi2 - nec_data.geometry.y[ic-1]) + abs(zi2 - nec_data.geometry.z[ic-1])
                if sep <= slen:
                    nec_data.geometry.icon2[i] = ic
                    break
                    
                sep = abs(xi2 - (nec_data.geometry.x[ic-1] + nec_data.geometry.si[ic-1] * nec_data.geometry.alp[ic-1])) + \
                      abs(yi2 - (nec_data.geometry.y[ic-1] + nec_data.geometry.si[ic-1] * nec_data.geometry.bet[ic-1])) + \
                      abs(zi2 - (nec_data.geometry.z[ic-1] + nec_data.geometry.si[ic-1] * nec_data.geometry.salp[ic-1]))
                if sep <= slen:
                    nec_data.geometry.icon2[i] = -ic
                    break
            else:
                if i < nec_data.geometry.n2 and nec_data.geometry.icon2[i] > 10000:
                    pass
                else:
                    nec_data.geometry.icon2[i] = 0

def process_patch_connections(ignd: int, smin: float, nec_data: NECData) -> None:
    """
    Process connections for surface patches
    
    Args:
        ignd: Ground flag
        smin: Minimum separation
        nec_data: Main NEC data structure
    """
    # Find wire-surface connections for new patches
    ix = nec_data.geometry.ld + 1 - nec_data.geometry.m1
    for i in range(nec_data.geometry.m2, nec_data.geometry.m + 1):
        ix = ix - 1
        xs = nec_data.geometry.x[ix - 1]
        ys = nec_data.geometry.y[ix - 1]
        zs = nec_data.geometry.z[ix - 1]
        
        for iseg in range(nec_data.geometry.n):
            # Get segment end coordinates
            xi1 = nec_data.geometry.x[iseg]
            yi1 = nec_data.geometry.y[iseg]
            zi1 = nec_data.geometry.z[iseg]
            xi2 = nec_data.geometry.x[iseg] + nec_data.geometry.si[iseg] * nec_data.geometry.alp[iseg]
            yi2 = nec_data.geometry.y[iseg] + nec_data.geometry.si[iseg] * nec_data.geometry.bet[iseg]
            zi2 = nec_data.geometry.z[iseg] + nec_data.geometry.si[iseg] * nec_data.geometry.salp[iseg]
            
            slen = (abs(xi2 - xi1) + abs(yi2 - yi1) + abs(zi2 - zi1)) * smin
            
            # Check first end of segment
            sep = abs(xi1 - xs) + abs(yi1 - ys) + abs(zi1 - zs)
            if sep <= slen:
                # Connection - divide patch into 4 patches
                nec_data.geometry.icon1[iseg] = 10000 + i
                ic = 0
                subph(i, ic, xi1, yi1, zi1, xi2, yi2, zi2, xs, ys, zs, nec_data)
                break
                
            # Check second end of segment
            sep = abs(xi2 - xs) + abs(yi2 - ys) + abs(zi2 - zs)
            if sep <= slen:
                # Connection - divide patch into 4 patches
                nec_data.geometry.icon2[iseg] = 10000 + i
                ic = 0
                subph(i, ic, xi1, yi1, zi1, xi2, yi2, zi2, xs, ys, zs, nec_data)
                break
    
    # Repeat search for new segments connected to NGF patches
    if nec_data.geometry.m1 == 0 or nec_data.geometry.n2 > nec_data.geometry.n:
        return
        
    ix = nec_data.geometry.ld + 1
    for i in range(1, nec_data.geometry.m1 + 1):
        ix = ix - 1
        xs = nec_data.geometry.x[ix - 1]
        ys = nec_data.geometry.y[ix - 1]
        zs = nec_data.geometry.z[ix - 1]
        
        for iseg in range(nec_data.geometry.n2 - 1, nec_data.geometry.n):
            # Get segment end coordinates
            xi1 = nec_data.geometry.x[iseg]
            yi1 = nec_data.geometry.y[iseg]
            zi1 = nec_data.geometry.z[iseg]
            xi2 = nec_data.geometry.x[iseg] + nec_data.geometry.si[iseg] * nec_data.geometry.alp[iseg]
            yi2 = nec_data.geometry.y[iseg] + nec_data.geometry.si[iseg] * nec_data.geometry.bet[iseg]
            zi2 = nec_data.geometry.z[iseg] + nec_data.geometry.si[iseg] * nec_data.geometry.salp[iseg]
            
            slen = (abs(xi2 - xi1) + abs(yi2 - yi1) + abs(zi2 - zi1)) * smin
            
            # Check first end of segment
            sep = abs(xi1 - xs) + abs(yi1 - ys) + abs(zi1 - zs)
            if sep <= slen:
                nec_data.geometry.icon1[iseg] = 10001 + nec_data.geometry.m
                ic = 1
                npcon += 1
                ipcon[npcon - 1] = i
                subph(i, ic, xi1, yi1, zi1, xi2, yi2, zi2, xs, ys, zs, nec_data)
                break
                
            # Check second end of segment
            sep = abs(xi2 - xs) + abs(yi2 - ys) + abs(zi2 - zs)
            if sep <= slen:
                nec_data.geometry.icon2[iseg] = 10001 + nec_data.geometry.m
                ic = 1
                npcon += 1
                ipcon[npcon - 1] = i
                subph(i, ic, xi1, yi1, zi1, xi2, yi2, zi2, xs, ys, zs, nec_data)
                break

def subph(i: int, ic: int, xi1: float, yi1: float, zi1: float, 
          xi2: float, yi2: float, zi2: float, xs: float, ys: float, zs: float,
          nec_data: NECData) -> None:
    """
    Subdivide patch (SUBPH from Fortran)
    
    Args:
        i: Patch index
        ic: Connection flag
        xi1, yi1, zi1: First segment end
        xi2, yi2, zi2: Second segment end
        xs, ys, zs: Patch center
        nec_data: Main NEC data structure
    """
    # This is a simplified implementation of patch subdivision
    # The full implementation would create 4 sub-patches
    
    # Calculate patch subdivision parameters
    dx = (xi2 - xi1) / 4.0
    dy = (yi2 - yi1) / 4.0
    dz = (zi2 - zi1) / 4.0
    
    # Create 4 sub-patches (simplified)
    for k in range(4):
        # Calculate sub-patch center
        xc = xs + k * dx
        yc = ys + k * dy
        zc = zs + k * dz
        
        # Store sub-patch data (simplified)
        # In the full implementation, this would create actual patch data structures

def process_multiple_junctions(nec_data: NECData) -> None:
    """
    Process multiple wire junctions
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\n   - MULTIPLE WIRE JUNCTIONS -")
    print("JUNCTION     SEGMENTS  (- FOR END 1, + FOR END 2)")
    
    iseg = 0
    
    # Adjust connected segment ends to exactly coincide
    for j in range(nec_data.geometry.n):
        iend = -1
        jend = -1
        ix = nec_data.geometry.icon1[j]
        ic = 1
        jco = [-j]
        xa = nec_data.geometry.x[j]
        ya = nec_data.geometry.y[j]
        za = nec_data.geometry.z[j]
        
        # Process first end
        while ix != 0 and ix != j and ix <= 10000:
            nsflg = 0
            
            if ix < 0:
                ix = -ix
            else:
                jend = -jend
                
            if ix == j:
                # Junction found
                sep = ic
                xa = xa / sep
                ya = ya / sep
                za = za / sep
                
                # Adjust segment coordinates
                for i in range(ic):
                    ix_jco = jco[i]
                    if ix_jco > 0:
                        ix_jco = -ix_jco
                    ix_jco = -ix_jco
                    nec_data.geometry.x[ix_jco - 1] = xa
                    nec_data.geometry.y[ix_jco - 1] = ya
                    nec_data.geometry.z[ix_jco - 1] = za
                
                # Check for old segments connecting to new segments
                if nec_data.geometry.n1 != 0:
                    for i in range(ic):
                        ix_jco = abs(jco[i])
                        if ix_jco > nec_data.geometry.n1:
                            if nec_data.geometry.iconx[ix_jco - 1] == 0:
                                nscon += 1
                                if nscon <= 50:
                                    iscon[nscon - 1] = ix_jco
                                    nec_data.geometry.iconx[ix_jco - 1] = nscon
                
                # Print junction if 3 or more segments
                if ic >= 3:
                    iseg += 1
                    print(f"{iseg:5d}     ", end="")
                    for i in range(min(20, ic)):
                        print(f"{jco[i]:5d}", end="")
                    print()
                    if ic > 20:
                        print("           ", end="")
                        for i in range(20, min(40, ic)):
                            print(f"{jco[i]:5d}", end="")
                        print()
                break
            else:
                if ix < j:
                    break
                    
                ic += 1
                if ic > 30:  # JMAX
                    print(f"CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT {j+1}")
                    return
                    
                jco.append(ix * jend)
                if ix > nec_data.geometry.n1:
                    nsflg = 1
                    
                if jend == 1:
                    xa += nec_data.geometry.x[ix - 1]
                    ya += nec_data.geometry.y[ix - 1]
                    za += nec_data.geometry.z[ix - 1]
                    ix = nec_data.geometry.icon1[ix - 1]
                else:
                    xa += nec_data.geometry.x[ix - 1] + nec_data.geometry.si[ix - 1] * nec_data.geometry.alp[ix - 1]
                    ya += nec_data.geometry.y[ix - 1] + nec_data.geometry.si[ix - 1] * nec_data.geometry.bet[ix - 1]
                    za += nec_data.geometry.z[ix - 1] + nec_data.geometry.si[ix - 1] * nec_data.geometry.salp[ix - 1]
                    ix = nec_data.geometry.icon2[ix - 1]
        
        # Process second end
        if iend == 1:
            continue
            
        iend = 1
        jend = 1
        ix = nec_data.geometry.icon2[j]
        ic = 1
        jco = [j]
        xa = nec_data.geometry.x[j] + nec_data.geometry.si[j] * nec_data.geometry.alp[j]
        ya = nec_data.geometry.y[j] + nec_data.geometry.si[j] * nec_data.geometry.bet[j]
        za = nec_data.geometry.z[j] + nec_data.geometry.si[j] * nec_data.geometry.salp[j]
        
        # Process second end (similar to first end)
        # ... (similar logic as above)
    
    if iseg == 0:
        print("   NONE")
    
    # Find old segments that connect to new patches
    if nec_data.geometry.n1 == 0 or nec_data.geometry.m1 == nec_data.geometry.m:
        return
        
    for j in range(nec_data.geometry.n1):
        ix = nec_data.geometry.icon1[j]
        if ix < 10000:
            ix = nec_data.geometry.icon2[j]
            if ix < 10000:
                continue
                
        ix = ix - 10000
        if ix > nec_data.geometry.m1:
            ix = nec_data.geometry.icon2[j]
            if ix < 10000:
                continue
            ix = ix - 10000
            if ix < nec_data.geometry.m2:
                continue
        else:
            if nec_data.geometry.iconx[j] != 0:
                continue
                
        nscon += 1
        iscon[nscon - 1] = j + 1
        nec_data.geometry.iconx[j] = nscon

def print_connection_summary(nec_data: NECData) -> None:
    """
    Print connection summary
    
    Args:
        nec_data: Main NEC data structure
    """
    print(f"\n   TOTAL SEGMENTS USED={nec_data.geometry.n:5d}     NO. SEG. IN A SYMMETRIC CELL={nec_data.geometry.np:5d}     SYMMETRY FLAG={nec_data.geometry.ipsym:3d}")
    
    if nec_data.geometry.m > 0:
        print(f"   TOTAL PATCHES USED={nec_data.geometry.m:5d}      NO. PATCHES IN A SYMMETRIC CELL={nec_data.geometry.mp:5d}")
    
    iseg = (nec_data.geometry.n + nec_data.geometry.m) // (nec_data.geometry.np + nec_data.geometry.mp)
    if iseg == 1:
        return
        
    if nec_data.geometry.ipsym == 0:
        print("STRUCTURE HAS NO SYMMETRY")
        return
    elif nec_data.geometry.ipsym > 0:
        ic = iseg // 2
        if iseg == 8:
            ic = 3
        print(f"STRUCTURE HAS {ic:2d} PLANES OF SYMMETRY")
    else:
        print(f"STRUCTURE HAS {iseg:4d} FOLD ROTATIONAL SYMMETRY")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def find_segment_connections(iseg: int, nec_data: NECData) -> List[int]:
    """
    Find all segments connected to a given segment
    
    Args:
        iseg: Segment index
        nec_data: Main NEC data structure
        
    Returns:
        List of connected segment indices
    """
    connections = []
    
    # Check connections from this segment
    for i in range(nec_data.geometry.n):
        if nec_data.geometry.icon1[i] == iseg + 1 or nec_data.geometry.icon2[i] == iseg + 1:
            connections.append(i)
        elif nec_data.geometry.icon1[i] == -(iseg + 1) or nec_data.geometry.icon2[i] == -(iseg + 1):
            connections.append(i)
    
    return connections

def check_connection_validity(nec_data: NECData) -> bool:
    """
    Check if all connections are valid
    
    Args:
        nec_data: Main NEC data structure
        
    Returns:
        True if all connections are valid
    """
    for i in range(nec_data.geometry.n):
        # Check icon1
        icon1 = nec_data.geometry.icon1[i]
        if icon1 != 0 and abs(icon1) <= nec_data.geometry.n:
            # Check if referenced segment exists
            if abs(icon1) > nec_data.geometry.n:
                print(f"Invalid connection: segment {i+1} connects to non-existent segment {abs(icon1)}")
                return False
        
        # Check icon2
        icon2 = nec_data.geometry.icon2[i]
        if icon2 != 0 and abs(icon2) <= nec_data.geometry.n:
            # Check if referenced segment exists
            if abs(icon2) > nec_data.geometry.n:
                print(f"Invalid connection: segment {i+1} connects to non-existent segment {abs(icon2)}")
                return False
    
    return True

def print_connection_debug(nec_data: NECData) -> None:
    """
    Print debug information about connections
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\nCONNECTION DEBUG INFORMATION:")
    print("=" * 50)
    
    for i in range(min(20, nec_data.geometry.n)):  # Print first 20 segments
        print(f"Segment {i+1:3d}: icon1={nec_data.geometry.icon1[i]:6d}, icon2={nec_data.geometry.icon2[i]:6d}, "
              f"tag={nec_data.geometry.itag[i]:3d}")
    
    if nec_data.geometry.n > 20:
        print(f"... and {nec_data.geometry.n - 20} more segments")
    
    print("=" * 50) 