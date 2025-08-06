"""
Data Generation Module for Python NEC-2 Implementation

This module contains the data generation functions converted from the Fortran code,
including DATAGN and related subroutines for input of geometry data.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from .constants import *
from .data_structures import *
from .geometry import *
from .input_output import *

# ============================================================================
# MAIN DATA GENERATION ROUTINE (DATAGN)
# ============================================================================

def datagn(nec_data: NECData) -> None:
    """
    Main routine for input of geometry data (DATAGN from Fortran)
    
    Args:
        nec_data: Main NEC data structure
    """
    # Initialize geometry parameters
    nec_data.geometry.ipsym = 0
    nwire = 0
    nec_data.geometry.n = 0
    nec_data.geometry.np = 0
    nec_data.geometry.m = 0
    nec_data.geometry.mp = 0
    nec_data.geometry.n1 = 0
    nec_data.geometry.n2 = 1
    nec_data.geometry.m1 = 0
    nec_data.geometry.m2 = 1
    isct = 0
    iphd = 0
    
    print("\n" + "="*60)
    print("STRUCTURE SPECIFICATION")
    print("="*60)
    print("COORDINATES MUST BE INPUT IN")
    print("METERS OR BE SCALED TO METERS")
    print("BEFORE STRUCTURE INPUT IS ENDED")
    print()
    
    print("WIRE                                                           NO. OF  FIRST  LAST   TAG")
    print("NO.        X1        Y1        Z1         X2        Y2        Z2      RADIUS  SEG.   SEG.   SEG.   NO.")
    
    # Main geometry input loop
    while True:
        try:
            # Read geometry data card
            card_data = read_geometry_card()
            if card_data is None:
                break
                
            gm, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad = card_data
            
            # Check if we have room for more geometry
            if nec_data.geometry.n + nec_data.geometry.m > nec_data.geometry.ld:
                print("NUMBER OF WIRE SEGMENTS AND SURFACE PATCHES EXCEEDS DIMENSION LIMIT.")
                return
            
            # Process different geometry card types
            if gm == 'GF':
                # Read numerical Green's function tape
                process_green_function_card(itg, nec_data)
            elif gm == 'GE':
                # End geometry input
                process_end_geometry_card(ns, nec_data)
                break
            elif gm == 'GW':
                # Generate wire
                nwire = process_wire_card(nwire, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            elif gm == 'GA':
                # Generate wire arc
                nwire = process_wire_arc_card(nwire, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            elif gm == 'GH':
                # Generate helix
                nwire = process_helix_card(nwire, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            elif gm == 'SP':
                # Generate surface patch
                process_surface_patch_card(itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            elif gm == 'GX':
                # Reflect structure
                process_reflect_card(itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            elif gm == 'GR':
                # Rotate structure
                process_rotate_card(itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            elif gm == 'GS':
                # Scale structure
                process_scale_card(xw1, nec_data)
            elif gm == 'GM':
                # Move structure
                process_move_card(itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, nec_data)
            else:
                print(f"GEOMETRY DATA CARD ERROR: {gm}")
                print(f"Card data: {card_data}")
                return
                
        except Exception as e:
            print(f"Error processing geometry card: {e}")
            break

def read_geometry_card() -> Optional[Tuple[str, int, int, float, float, float, float, float, float, float]]:
    """
    Read a geometry card from input
    
    Returns:
        Tuple of (card_type, tag, segments, x1, y1, z1, x2, y2, z2, radius) or None if EOF
    """
    # This is a simplified implementation
    # In the full implementation, this would read from the actual input file
    
    # Placeholder - would read actual geometry card
    # For now, return None to end geometry input
    return None

def process_wire_card(nwire: int, itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                     xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> int:
    """
    Process wire card (GW from Fortran)
    
    Args:
        nwire: Wire number
        itg: Tag number
        ns: Number of segments
        xw1, yw1, zw1: First end coordinates
        xw2, yw2, zw2: Second end coordinates
        rad: Wire radius
        nec_data: Main NEC data structure
        
    Returns:
        Updated wire number
    """
    nwire += 1
    i1 = nec_data.geometry.n + 1
    i2 = nec_data.geometry.n + ns
    
    print(f"{nwire:5d} {xw1:11.5f} {yw1:11.5f} {zw1:11.5f} {xw2:11.5f} {yw2:11.5f} {zw2:11.5f} "
          f"{rad:11.5f} {ns:5d} {i1:5d} {i2:5d} {itg:5d}")
    
    if rad == 0.0:
        # Tapered wire - read additional card
        taper_data = read_taper_card()
        if taper_data is None:
            print("GEOMETRY DATA CARD ERROR")
            return nwire
            
        xs1, ys1, zs1, dummy, dummy, dummy, dummy = taper_data
        print(f"ABOVE WIRE IS TAPERED.  SEG. LENGTH RATIO = {xs1:.5f}")
        print(f"RADIUS FROM {ys1:.5f} TO {zs1:.5f}")
        rad = ys1
        ys1 = (zs1 / ys1) ** (1.0 / (ns - 1.0))
    else:
        xs1 = 1.0
        ys1 = 1.0
    
    # Generate wire segments
    generate_wire_segments(nec_data, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, xs1, ys1, itg)
    
    return nwire

def process_wire_arc_card(nwire: int, itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                         xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> int:
    """
    Process wire arc card (GA from Fortran)
    
    Args:
        nwire: Wire number
        itg: Tag number
        ns: Number of segments
        xw1, yw1, zw1: Center coordinates
        xw2, yw2, zw2: Arc parameters
        rad: Wire radius
        nec_data: Main NEC data structure
        
    Returns:
        Updated wire number
    """
    nwire += 1
    i1 = nec_data.geometry.n + 1
    i2 = nec_data.geometry.n + ns
    
    print(f"{nwire:5d} ARC RADIUS = {xw1:.5f} FROM {yw1:.3f} TO {zw1:.3f} DEGREES "
          f"{xw2:11.5f} {ns:5d} {i1:5d} {i2:5d} {itg:5d}")
    
    # Generate arc segments
    generate_arc_segments(nec_data, ns, xw1, yw1, zw1, xw2, itg)
    
    return nwire

def process_helix_card(nwire: int, itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                      xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> int:
    """
    Process helix card (GH from Fortran)
    
    Args:
        nwire: Wire number
        itg: Tag number
        ns: Number of segments
        xw1, yw1, zw1: Center coordinates
        xw2, yw2, zw2: Helix parameters
        rad: Wire radius
        nec_data: Main NEC data structure
        
    Returns:
        Updated wire number
    """
    nwire += 1
    i1 = nec_data.geometry.n + 1
    i2 = nec_data.geometry.n + ns
    
    print(f"HELIX STRUCTURE-   AXIAL SPACING BETWEEN TURNS = {xw1:.3f} TOTAL AXIAL LENGTH = {yw1:.3f}")
    print(f"{nwire:5d} RADIUS OF HELIX = {zw1:.3f} {xw2:.3f} {yw2:.3f} {zw2:.3f} {rad:11.5f} {ns:8d} {i1:5d} {i2:5d} {itg:5d}")
    
    # Generate helix segments
    generate_helix_segments(nec_data, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, itg)
    
    return nwire

def process_surface_patch_card(itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                              xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> None:
    """
    Process surface patch card (SP from Fortran)
    
    Args:
        itg: Tag number
        ns: Number of segments
        xw1, yw1, zw1: First corner coordinates
        xw2, yw2, zw2: Second corner coordinates
        rad: Patch area
        nec_data: Main NEC data structure
    """
    i1 = nec_data.geometry.m + 1
    ns = ns + 1
    
    if itg != 0:
        print(f"{i1:5d}P {xw1:10.5f} {yw1:11.5f} {zw1:11.5f} {xw2:11.5f} {yw2:11.5f} {zw2:11.5f}")
        if ns == 2 or ns == 4:
            isct = 1
        if ns > 1:
            # Read additional corner coordinates
            corner_data = read_corner_card()
            if corner_data is not None:
                x3, y3, z3, x4, y4, z4, dummy = corner_data
                print(f"      {x3:11.5f} {y3:11.5f} {z3:11.5f} {x4:11.5f} {y4:11.5f} {z4:11.5f}")
                generate_surface_patches(nec_data, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4, itg)
    else:
        print("PATCH DATA ERROR")

def process_reflect_card(itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                        xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> None:
    """
    Process reflect card (GX from Fortran)
    
    Args:
        itg: Tag number
        ns: Symmetry parameters
        xw1, yw1, zw1: Unused
        xw2, yw2, zw2: Unused
        rad: Unused
        nec_data: Main NEC data structure
    """
    iy = ns // 10
    iz = ns - iy * 10
    ix = iy // 10
    iy = iy - ix * 10
    
    if ix != 0:
        ix = 1
    if iy != 0:
        iy = 1
    if iz != 0:
        iz = 1
    
    print(f"STRUCTURE REFLECTED ALONG THE AXES {'X' if ix else ' '}{'Y' if iy else ' '}{'Z' if iz else ' '}.  TAGS INCREMENTED BY {itg}")
    
    # Reflect structure
    reflect_structure(ix, iy, iz, itg, ns, nec_data)

def process_rotate_card(itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                       xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> None:
    """
    Process rotate card (GR from Fortran)
    
    Args:
        itg: Tag number
        ns: Number of rotations
        xw1, yw1, zw1: Unused
        xw2, yw2, zw2: Unused
        rad: Unused
        nec_data: Main NEC data structure
    """
    print(f"STRUCTURE ROTATED ABOUT Z-AXIS {ns} TIMES.  LABELS INCREMENTED BY {itg}")
    
    # Rotate structure
    rotate_structure(ns, itg, nec_data)

def process_scale_card(xw1: float, nec_data: NECData) -> None:
    """
    Process scale card (GS from Fortran)
    
    Args:
        xw1: Scale factor
        nec_data: Main NEC data structure
    """
    # Scale wire segments
    if nec_data.geometry.n >= nec_data.geometry.n2:
        for i in range(nec_data.geometry.n2 - 1, nec_data.geometry.n):
            nec_data.geometry.x[i] *= xw1
            nec_data.geometry.y[i] *= xw1
            nec_data.geometry.z[i] *= xw1
            nec_data.geometry.si[i] *= xw1
            nec_data.geometry.alp[i] *= xw1
            nec_data.geometry.bet[i] *= xw1
            nec_data.geometry.bi[i] *= xw1
    
    # Scale surface patches
    if nec_data.geometry.m >= nec_data.geometry.m2:
        yw1 = xw1 * xw1
        ix = nec_data.geometry.ld + 1 - nec_data.geometry.m
        iy = nec_data.geometry.ld - nec_data.geometry.m1
        for i in range(ix - 1, iy):
            nec_data.geometry.x[i] *= xw1
            nec_data.geometry.y[i] *= xw1
            nec_data.geometry.z[i] *= xw1
            nec_data.geometry.bi[i] *= yw1
    
    print(f"STRUCTURE SCALED BY FACTOR {xw1:.5f}")

def process_move_card(itg: int, ns: int, xw1: float, yw1: float, zw1: float,
                     xw2: float, yw2: float, zw2: float, rad: float, nec_data: NECData) -> None:
    """
    Process move card (GM from Fortran)
    
    Args:
        itg: Tag number
        ns: Number of copies
        xw1, yw1, zw1: Translation parameters
        xw2, yw2, zw2: Rotation parameters
        rad: Number of copies
        nec_data: Main NEC data structure
    """
    print(f"THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS -")
    print(f"{itg:6d} {ns:5d} {xw1:10.5f} {yw1:10.5f} {zw1:10.5f} {xw2:10.5f} {yw2:10.5f} {zw2:10.5f} {rad:10.5f}")
    
    # Convert angles to radians
    xw1 *= TA
    yw1 *= TA
    zw1 *= TA
    
    # Move structure
    move_structure(xw1, yw1, zw1, xw2, yw2, zw2, int(rad + 0.5), ns, itg, nec_data)

def process_green_function_card(itg: int, nec_data: NECData) -> None:
    """
    Process Green's function card (GF from Fortran)
    
    Args:
        itg: Tag number
        nec_data: Main NEC data structure
    """
    if nec_data.geometry.n + nec_data.geometry.m != 0:
        print("ERROR - GF MUST BE FIRST GEOMETRY DATA CARD")
        return
    
    # Read Green's function file
    read_green_function_file(itg, nec_data)
    
    # Save parameters
    nec_data.geometry.np = nec_data.geometry.n
    nec_data.geometry.mp = nec_data.geometry.m
    nec_data.geometry.ipsym = 0

def process_end_geometry_card(ns: int, nec_data: NECData) -> None:
    """
    Process end geometry card (GE from Fortran)
    
    Args:
        ns: Plot flag
        nec_data: Main NEC data structure
    """
    if ns != 0:
        nec_data.iplp1 = 1
        nec_data.iplp2 = 1
    
    ix = nec_data.geometry.n1 + nec_data.geometry.m1
    if ix == 0:
        # Set up connections
        setup_connections(nec_data)
    else:
        # Restore parameters for NGF
        nec_data.geometry.np = nec_data.geometry.n
        nec_data.geometry.mp = nec_data.geometry.m
        nec_data.geometry.ipsym = 0
        setup_connections(nec_data)
        # Restore NGF parameters
        # (This would restore the saved parameters)
    
    if nec_data.geometry.n + nec_data.geometry.m > nec_data.geometry.ld:
        print("NUMBER OF WIRE SEGMENTS AND SURFACE PATCHES EXCEEDS DIMENSION LIMIT.")
        return
    
    # Print geometry summary
    print_geometry_summary(nec_data)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def read_taper_card() -> Optional[Tuple[float, float, float, float, float, float, float]]:
    """
    Read taper card for wire
    
    Returns:
        Taper parameters or None if error
    """
    # This is a simplified implementation
    # In the full implementation, this would read from the actual input file
    return None

def read_corner_card() -> Optional[Tuple[float, float, float, float, float, float, float]]:
    """
    Read corner coordinates for surface patch
    
    Returns:
        Corner coordinates or None if error
    """
    # This is a simplified implementation
    # In the full implementation, this would read from the actual input file
    return None

def reflect_structure(ix: int, iy: int, iz: int, itg: int, ns: int, nec_data: NECData) -> None:
    """
    Reflect structure along specified axes
    
    Args:
        ix: X-axis reflection flag
        iy: Y-axis reflection flag
        iz: Z-axis reflection flag
        itg: Tag increment
        ns: Symmetry parameters
        nec_data: Main NEC data structure
    """
    # This is a simplified implementation
    # The full implementation would create reflected copies of segments and patches
    
    # Update symmetry flag
    if ix != 0 or iy != 0 or iz != 0:
        nec_data.geometry.ipsym = 1
    
    # Reflect wire segments
    if iz != 0 and nec_data.geometry.n >= nec_data.geometry.n2:
        # Reflect along Z-axis
        for i in range(nec_data.geometry.n2 - 1, nec_data.geometry.n):
            # Create reflected segment
            # (Simplified - would create actual reflected segments)
            pass
    
    # Reflect surface patches
    if iz != 0 and nec_data.geometry.m >= nec_data.geometry.m2:
        # Reflect patches along Z-axis
        # (Simplified - would create actual reflected patches)
        pass

def rotate_structure(ns: int, itg: int, nec_data: NECData) -> None:
    """
    Rotate structure to form cylindrical structure
    
    Args:
        ns: Number of rotations
        itg: Tag increment
        nec_data: Main NEC data structure
    """
    # This is a simplified implementation
    # The full implementation would create rotated copies of segments and patches
    
    # Update symmetry flag
    nec_data.geometry.ipsym = -1
    
    # Calculate rotation angle
    sam = 2 * PI / ns
    cs = np.cos(sam)
    ss = np.sin(sam)
    
    # Rotate wire segments
    if nec_data.geometry.n >= nec_data.geometry.n2:
        # Create rotated segments
        # (Simplified - would create actual rotated segments)
        pass
    
    # Rotate surface patches
    if nec_data.geometry.m >= nec_data.geometry.m2:
        # Create rotated patches
        # (Simplified - would create actual rotated patches)
        pass

def move_structure(xw1: float, yw1: float, zw1: float, xw2: float, yw2: float, zw2: float,
                  ncopies: int, ns: int, itg: int, nec_data: NECData) -> None:
    """
    Move structure or reproduce in new positions
    
    Args:
        xw1, yw1, zw1: Translation parameters
        xw2, yw2, zw2: Rotation parameters
        ncopies: Number of copies
        ns: Number of copies
        itg: Tag increment
        nec_data: Main NEC data structure
    """
    # This is a simplified implementation
    # The full implementation would create moved copies of segments and patches
    
    for i in range(ncopies):
        # Calculate new position
        # (Simplified - would calculate actual new positions)
        pass
        
        # Create moved segments and patches
        # (Simplified - would create actual moved geometry)
        pass

def read_green_function_file(itg: int, nec_data: NECData) -> None:
    """
    Read numerical Green's function file
    
    Args:
        itg: Tag number
        nec_data: Main NEC data structure
    """
    # This is a simplified implementation
    # The full implementation would read the actual NGF file
    
    print("READING NUMERICAL GREEN'S FUNCTION FILE")
    print("(This is a placeholder - actual file reading not implemented)")

def print_geometry_summary(nec_data: NECData) -> None:
    """
    Print geometry summary
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\n" + "="*60)
    print("GEOMETRY SUMMARY")
    print("="*60)
    print(f"Total wire segments: {nec_data.geometry.n}")
    print(f"Total surface patches: {nec_data.geometry.m}")
    print(f"Wire segments in symmetric cell: {nec_data.geometry.np}")
    print(f"Surface patches in symmetric cell: {nec_data.geometry.mp}")
    print(f"Symmetry flag: {nec_data.geometry.ipsym}")
    
    if nec_data.geometry.n > 0:
        print("\nWIRE SEGMENT DATA:")
        print("SEG.  COORDINATES OF SEG. CENTER  SEG.  ORIENTATION ANGLES  WIRE  CONNECTION DATA  TAG")
        print("NO.       X        Y        Z   LENGTH    ALPHA     BETA  RADIUS    I-   I   I+   NO.")
        
        for i in range(min(20, nec_data.geometry.n)):  # Print first 20 segments
            xc = (nec_data.geometry.x[i] + nec_data.geometry.x[i] + nec_data.geometry.si[i] * nec_data.geometry.alp[i]) / 2
            yc = (nec_data.geometry.y[i] + nec_data.geometry.y[i] + nec_data.geometry.si[i] * nec_data.geometry.bet[i]) / 2
            zc = (nec_data.geometry.z[i] + nec_data.geometry.z[i] + nec_data.geometry.si[i] * nec_data.geometry.salp[i]) / 2
            
            print(f"{i+1:4d} {xc:9.5f} {yc:9.5f} {zc:9.5f} {nec_data.geometry.si[i]:9.5f} "
                  f"{nec_data.geometry.alp[i]:9.5f} {nec_data.geometry.bet[i]:9.5f} {nec_data.geometry.bi[i]:9.5f} "
                  f"{nec_data.geometry.icon1[i]:5d} {i+1:3d} {nec_data.geometry.icon2[i]:3d} {nec_data.geometry.itag[i]:5d}")
        
        if nec_data.geometry.n > 20:
            print(f"... and {nec_data.geometry.n - 20} more segments")
    
    if nec_data.geometry.m > 0:
        print("\nSURFACE PATCH DATA:")
        print("PATCH  COORD. OF PATCH CENTER  UNIT NORMAL VECTOR  PATCH  COMPONENTS OF UNIT TANGENT VECTORS")
        print("NO.        X        Y        Z       X      Y      Z   AREA     X1     Y1     Z1     X2     Y2     Z2")
        
        for i in range(min(10, nec_data.geometry.m)):  # Print first 10 patches
            j = nec_data.geometry.ld - i
            print(f"{i+1:4d} {nec_data.geometry.x[j]:9.5f} {nec_data.geometry.y[j]:9.5f} {nec_data.geometry.z[j]:9.5f} "
                  f"{nec_data.geometry.alp[j]:8.4f} {nec_data.geometry.bet[j]:8.4f} {nec_data.geometry.salp[j]:8.4f} "
                  f"{nec_data.geometry.bi[j]:9.5f} {nec_data.geometry.si[j]:8.4f} {nec_data.geometry.alp[j]:8.4f} "
                  f"{nec_data.geometry.bet[j]:8.4f} {nec_data.geometry.salp[j]:8.4f} {nec_data.geometry.si[j]:8.4f} "
                  f"{nec_data.geometry.alp[j]:8.4f} {nec_data.geometry.bet[j]:8.4f}")
        
        if nec_data.geometry.m > 10:
            print(f"... and {nec_data.geometry.m - 10} more patches")
    
    print("="*60) 