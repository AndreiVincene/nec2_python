"""
Geometry Processing Module for Python NEC-2 Implementation

This module contains all geometry-related functions converted from the Fortran code,
including DATAGN, WIRE, ARC, HELIX, and other geometry processing subroutines.
"""

import numpy as np
from typing import Tuple, List, Optional
from .constants import *
from .data_structures import *

# ============================================================================
# MAIN GEOMETRY DATA GENERATION (DATAGN)
# ============================================================================

def datagn(nec_data: NECData) -> None:
    """
    Main geometry data generation routine (DATAGN from Fortran)
    
    This is the main routine for input of geometry data, corresponding to
    the DATAGN subroutine in the original Fortran code.
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
    
    # Read geometry data cards and process them
    while True:
        try:
            # Read geometry card
            card_data = read_geometry_card()
            if card_data is None:
                break
                
            gm, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad = card_data
            
            # Check if we exceed maximum segments
            if nec_data.geometry.n + nec_data.geometry.m > nec_data.geometry.ld:
                raise ValueError("Number of segments exceeds maximum allowed")
            
            # Process different geometry card types
            if gm == 'GF':  # Surface file
                process_surface_file(nec_data)
                continue
                
            if iphd == 1:
                pass  # Already printed header
            else:
                print_geometry_header()
                iphd = 1
            
            # Process geometry based on card type
            if gm == 'SC':  # Surface cone
                isct = 0
            elif gm == 'GW':  # Wire
                nwire = process_wire(nec_data, nwire, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GX':  # Wire arc
                nwire = process_wire_arc(nec_data, nwire, itg, ns, xw1, yw1, zw1, xw2)
            elif gm == 'GR':  # Wire grid
                process_wire_grid(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GS':  # Wire sphere
                process_wire_sphere(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'SP':  # Surface patch
                process_surface_patch(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'SM':  # Surface move
                process_surface_move(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GE':  # Wire ellipse
                process_wire_ellipse(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GM':  # Wire move
                process_wire_move(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GA':  # Surface arc
                process_surface_arc(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GC':  # Surface cylinder
                process_surface_cylinder(nec_data, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            elif gm == 'GH':  # Helix
                nwire = process_helix(nec_data, nwire, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad)
            else:
                raise ValueError(f"Unknown geometry card type: {gm}")
                
        except EOFError:
            break
    
    # Set up connection data
    setup_connections(nec_data)
    
    # Print geometry summary
    print_geometry_summary(nec_data)

def read_geometry_card() -> Optional[Tuple[str, int, int, float, float, float, float, float, float, float]]:
    """
    Read a geometry card from input (simplified version)
    
    Returns:
        Tuple of (card_type, tag, segments, x1, y1, z1, x2, y2, z2, radius)
        or None if end of input
    """
    # This is a simplified version - in the full implementation,
    # this would read from the actual input file
    return None

def print_geometry_header() -> None:
    """Print geometry section header"""
    print("=" * 60)
    print("GEOMETRY DATA")
    print("=" * 60)

def print_geometry_summary(nec_data: NECData) -> None:
    """Print geometry summary"""
    print(f"\nTotal wire segments: {nec_data.geometry.n}")
    print(f"Wire segments in symmetric cell: {nec_data.geometry.np}")
    print(f"Total surface patches: {nec_data.geometry.m}")
    print(f"Surface patches in symmetric cell: {nec_data.geometry.mp}")
    print(f"Symmetry flag: {nec_data.geometry.ipsym}")

# ============================================================================
# WIRE GEOMETRY PROCESSING
# ============================================================================

def process_wire(nec_data: NECData, nwire: int, itg: int, ns: int, 
                xw1: float, yw1: float, zw1: float, 
                xw2: float, yw2: float, zw2: float, rad: float) -> int:
    """
    Process wire geometry card (GW)
    
    Args:
        nec_data: NEC data structure
        nwire: Current wire count
        itg: Wire tag
        ns: Number of segments
        xw1, yw1, zw1: Start coordinates
        xw2, yw2, zw2: End coordinates
        rad: Wire radius
        
    Returns:
        Updated wire count
    """
    nwire += 1
    i1 = nec_data.geometry.n + 1
    i2 = nec_data.geometry.n + ns
    
    print(f"Wire {nwire}: ({xw1:.3f}, {yw1:.3f}, {zw1:.3f}) to ({xw2:.3f}, {yw2:.3f}, {zw2:.3f})")
    print(f"  Radius: {rad:.6f}, Segments: {ns}, Tag: {itg}, Segments {i1}-{i2}")
    
    # Generate wire segments
    if rad == 0.0:
        # Variable radius wire
        # Read additional card for radius specification
        pass
    else:
        # Constant radius wire
        generate_wire_segments(nec_data, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad, itg)
    
    return nwire

def generate_wire_segments(nec_data: NECData, ns: int, x1: float, y1: float, z1: float,
                          x2: float, y2: float, z2: float, rad: float, tag: int) -> None:
    """
    Generate wire segments for a straight wire
    
    Args:
        nec_data: NEC data structure
        ns: Number of segments
        x1, y1, z1: Start coordinates
        x2, y2, z2: End coordinates
        rad: Wire radius
        tag: Wire tag
    """
    # Calculate segment length
    dx = (x2 - x1) / ns
    dy = (y2 - y1) / ns
    dz = (z2 - z1) / ns
    
    # Calculate direction cosines
    length = np.sqrt(dx**2 + dy**2 + dz**2)
    if length > 0:
        cab = dx / length  # Alpha (cosine)
        sab = dy / length  # Beta (sine)
        salp = dz / length  # Gamma
    
    # Generate segments
    for i in range(ns):
        seg_idx = nec_data.geometry.n + i
        
        # Segment center coordinates
        nec_data.geometry.x[seg_idx] = x1 + (i + 0.5) * dx
        nec_data.geometry.y[seg_idx] = y1 + (i + 0.5) * dy
        nec_data.geometry.z[seg_idx] = z1 + (i + 0.5) * dz
        
        # Segment length
        nec_data.geometry.si[seg_idx] = length
        
        # Wire radius
        nec_data.geometry.bi[seg_idx] = rad
        
        # Direction cosines
        nec_data.geometry.alp[seg_idx] = cab
        nec_data.geometry.bet[seg_idx] = sab
        
        # Tag
        nec_data.geometry.itag[seg_idx] = tag
    
    # Update segment count
    nec_data.geometry.n += ns
    nec_data.geometry.np = nec_data.geometry.n

def process_wire_arc(nec_data: NECData, nwire: int, itg: int, ns: int,
                    xc: float, yc: float, zc: float, radius: float) -> int:
    """
    Process wire arc geometry card (GX)
    
    Args:
        nec_data: NEC data structure
        nwire: Current wire count
        itg: Wire tag
        ns: Number of segments
        xc, yc, zc: Center coordinates
        radius: Arc radius
        
    Returns:
        Updated wire count
    """
    nwire += 1
    i1 = nec_data.geometry.n + 1
    i2 = nec_data.geometry.n + ns
    
    print(f"Wire Arc {nwire}: Center ({xc:.3f}, {yc:.3f}, {zc:.3f}), Radius: {radius:.3f}")
    print(f"  Segments: {ns}, Tag: {itg}, Segments {i1}-{i2}")
    
    # Generate arc segments
    generate_arc_segments(nec_data, ns, xc, yc, zc, radius, itg)
    
    return nwire

def generate_arc_segments(nec_data: NECData, ns: int, xc: float, yc: float, zc: float,
                         radius: float, tag: int) -> None:
    """
    Generate wire segments for an arc
    
    Args:
        nec_data: NEC data structure
        ns: Number of segments
        xc, yc, zc: Center coordinates
        radius: Arc radius
        tag: Wire tag
    """
    # Arc angle increment
    dtheta = 2 * PI / ns
    
    for i in range(ns):
        seg_idx = nec_data.geometry.n + i
        
        # Calculate segment endpoints
        theta1 = i * dtheta
        theta2 = (i + 1) * dtheta
        
        x1 = xc + radius * np.cos(theta1)
        y1 = yc + radius * np.sin(theta1)
        z1 = zc
        
        x2 = xc + radius * np.cos(theta2)
        y2 = yc + radius * np.sin(theta2)
        z2 = zc
        
        # Segment center
        nec_data.geometry.x[seg_idx] = (x1 + x2) / 2
        nec_data.geometry.y[seg_idx] = (y1 + y2) / 2
        nec_data.geometry.z[seg_idx] = (z1 + z2) / 2
        
        # Segment length
        nec_data.geometry.si[seg_idx] = radius * dtheta
        
        # Wire radius (default)
        nec_data.geometry.bi[seg_idx] = radius * 0.001  # 0.1% of arc radius
        
        # Direction cosines
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        length = np.sqrt(dx**2 + dy**2 + dz**2)
        
        if length > 0:
            nec_data.geometry.alp[seg_idx] = dx / length
            nec_data.geometry.bet[seg_idx] = dy / length
        else:
            nec_data.geometry.alp[seg_idx] = 0.0
            nec_data.geometry.bet[seg_idx] = 0.0
        
        # Tag
        nec_data.geometry.itag[seg_idx] = tag
    
    # Update segment count
    nec_data.geometry.n += ns
    nec_data.geometry.np = nec_data.geometry.n

def process_helix(nec_data: NECData, nwire: int, itg: int, ns: int,
                 xc: float, yc: float, zc: float, 
                 xh: float, yh: float, zh: float, rad: float) -> int:
    """
    Process helix geometry card (GH)
    
    Args:
        nec_data: NEC data structure
        nwire: Current wire count
        itg: Wire tag
        ns: Number of segments
        xc, yc, zc: Center coordinates
        xh, yh, zh: Helix parameters
        rad: Wire radius
        
    Returns:
        Updated wire count
    """
    nwire += 1
    i1 = nec_data.geometry.n + 1
    i2 = nec_data.geometry.n + ns
    
    print(f"Helix {nwire}: Center ({xc:.3f}, {yc:.3f}, {zc:.3f})")
    print(f"  Parameters: ({xh:.3f}, {yh:.3f}, {zh:.3f}), Radius: {rad:.3f}")
    print(f"  Segments: {ns}, Tag: {itg}, Segments {i1}-{i2}")
    
    # Generate helix segments
    generate_helix_segments(nec_data, ns, xc, yc, zc, xh, yh, zh, rad, itg)
    
    return nwire

def generate_helix_segments(nec_data: NECData, ns: int, xc: float, yc: float, zc: float,
                           xh: float, yh: float, zh: float, rad: float, tag: int) -> None:
    """
    Generate wire segments for a helix
    
    Args:
        nec_data: NEC data structure
        ns: Number of segments
        xc, yc, zc: Center coordinates
        xh, yh, zh: Helix parameters
        rad: Wire radius
        tag: Wire tag
    """
    # Helix parameters
    radius = np.sqrt(xh**2 + yh**2)
    pitch = zh
    turns = 1.0  # Default to 1 turn
    
    # Angle increment
    dtheta = 2 * PI * turns / ns
    
    for i in range(ns):
        seg_idx = nec_data.geometry.n + i
        
        # Calculate segment endpoints
        theta1 = i * dtheta
        theta2 = (i + 1) * dtheta
        
        z1 = zc + (i * pitch / ns)
        z2 = zc + ((i + 1) * pitch / ns)
        
        x1 = xc + radius * np.cos(theta1)
        y1 = yc + radius * np.sin(theta1)
        
        x2 = xc + radius * np.cos(theta2)
        y2 = yc + radius * np.sin(theta2)
        
        # Segment center
        nec_data.geometry.x[seg_idx] = (x1 + x2) / 2
        nec_data.geometry.y[seg_idx] = (y1 + y2) / 2
        nec_data.geometry.z[seg_idx] = (z1 + z2) / 2
        
        # Segment length
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        nec_data.geometry.si[seg_idx] = np.sqrt(dx**2 + dy**2 + dz**2)
        
        # Wire radius
        nec_data.geometry.bi[seg_idx] = rad
        
        # Direction cosines
        length = nec_data.geometry.si[seg_idx]
        if length > 0:
            nec_data.geometry.alp[seg_idx] = dx / length
            nec_data.geometry.bet[seg_idx] = dy / length
        else:
            nec_data.geometry.alp[seg_idx] = 0.0
            nec_data.geometry.bet[seg_idx] = 0.0
        
        # Tag
        nec_data.geometry.itag[seg_idx] = tag
    
    # Update segment count
    nec_data.geometry.n += ns
    nec_data.geometry.np = nec_data.geometry.n

# ============================================================================
# SURFACE PATCH PROCESSING
# ============================================================================

def process_surface_patch(nec_data: NECData, itg: int, ns: int,
                         x1: float, y1: float, z1: float,
                         x2: float, y2: float, z2: float, rad: float) -> None:
    """
    Process surface patch geometry card (SP)
    
    Args:
        nec_data: NEC data structure
        itg: Patch tag
        ns: Number of segments
        x1, y1, z1: First corner coordinates
        x2, y2, z2: Second corner coordinates
        rad: Patch parameter
    """
    i1 = nec_data.geometry.m + 1
    i2 = nec_data.geometry.m + ns
    
    print(f"Surface Patch: ({x1:.3f}, {y1:.3f}, {z1:.3f}) to ({x2:.3f}, {y2:.3f}, {z2:.3f})")
    print(f"  Segments: {ns}, Tag: {itg}, Patches {i1}-{i2}")
    
    # Generate surface patches
    generate_surface_patches(nec_data, ns, x1, y1, z1, x2, y2, z2, itg)

def generate_surface_patches(nec_data: NECData, ns: int,
                            x1: float, y1: float, z1: float,
                            x2: float, y2: float, z2: float, tag: int) -> None:
    """
    Generate surface patches
    
    Args:
        nec_data: NEC data structure
        ns: Number of segments
        x1, y1, z1: First corner coordinates
        x2, y2, z2: Second corner coordinates
        tag: Patch tag
    """
    # Calculate patch dimensions
    dx = (x2 - x1) / ns
    dy = (y2 - y1) / ns
    dz = (z2 - z1) / ns
    
    for i in range(ns):
        patch_idx = nec_data.geometry.m + i
        
        # Patch center
        nec_data.geometry.x[patch_idx] = x1 + (i + 0.5) * dx
        nec_data.geometry.y[patch_idx] = y1 + (i + 0.5) * dy
        nec_data.geometry.z[patch_idx] = z1 + (i + 0.5) * dz
        
        # Patch area (simplified)
        nec_data.geometry.bi[patch_idx] = np.sqrt(dx**2 + dy**2 + dz**2)
        
        # Direction cosines (normal to patch)
        length = nec_data.geometry.bi[patch_idx]
        if length > 0:
            nec_data.geometry.alp[patch_idx] = dx / length
            nec_data.geometry.bet[patch_idx] = dy / length
        else:
            nec_data.geometry.alp[patch_idx] = 0.0
            nec_data.geometry.bet[patch_idx] = 0.0
        
        # Tag
        nec_data.geometry.itag[patch_idx] = tag
    
    # Update patch count
    nec_data.geometry.m += ns
    nec_data.geometry.mp = nec_data.geometry.m

# ============================================================================
# CONNECTION SETUP
# ============================================================================

def setup_connections(nec_data: NECData) -> None:
    """
    Set up segment connection data (CONECT from Fortran)
    
    This function sets up the connection data in arrays ICON1 and ICON2
    by searching for segment ends that are in contact.
    """
    # Initialize connection arrays
    nec_data.geometry.icon1.fill(0)
    nec_data.geometry.icon2.fill(0)
    nec_data.geometry.iconx.fill(0)
    
    # Connection tolerance
    smin = 1e-3
    
    # Find connections for wire segments
    for i in range(nec_data.geometry.n):
        # Segment endpoints
        x1 = nec_data.geometry.x[i] - 0.5 * nec_data.geometry.si[i] * nec_data.geometry.alp[i]
        y1 = nec_data.geometry.y[i] - 0.5 * nec_data.geometry.si[i] * nec_data.geometry.bet[i]
        z1 = nec_data.geometry.z[i] - 0.5 * nec_data.geometry.si[i] * np.sqrt(1 - nec_data.geometry.alp[i]**2 - nec_data.geometry.bet[i]**2)
        
        x2 = nec_data.geometry.x[i] + 0.5 * nec_data.geometry.si[i] * nec_data.geometry.alp[i]
        y2 = nec_data.geometry.y[i] + 0.5 * nec_data.geometry.si[i] * nec_data.geometry.bet[i]
        z2 = nec_data.geometry.z[i] + 0.5 * nec_data.geometry.si[i] * np.sqrt(1 - nec_data.geometry.alp[i]**2 - nec_data.geometry.bet[i]**2)
        
        # Find connections for end 1
        for j in range(nec_data.geometry.n):
            if i == j:
                continue
                
            # Check distance to other segment endpoints
            sep1 = abs(x1 - nec_data.geometry.x[j]) + abs(y1 - nec_data.geometry.y[j]) + abs(z1 - nec_data.geometry.z[j])
            sep2 = abs(x2 - nec_data.geometry.x[j]) + abs(y2 - nec_data.geometry.y[j]) + abs(z2 - nec_data.geometry.z[j])
            
            slen = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2) * smin
            
            if sep1 <= slen:
                nec_data.geometry.icon1[i] = j + 1
            if sep2 <= slen:
                nec_data.geometry.icon2[i] = j + 1
    
    print_connection_summary(nec_data)

def print_connection_summary(nec_data: NECData) -> None:
    """Print connection summary"""
    connections = 0
    for i in range(nec_data.geometry.n):
        if nec_data.geometry.icon1[i] != 0 or nec_data.geometry.icon2[i] != 0:
            connections += 1
    
    print(f"\nTotal segments: {nec_data.geometry.n}")
    print(f"Connected segments: {connections}")
    print(f"Symmetry flag: {nec_data.geometry.ipsym}")

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def atgn2(x: float, y: float) -> float:
    """
    Arctangent function modified to return 0 when x=y=0
    (ATGN2 from Fortran)
    """
    if x == 0.0 and y == 0.0:
        return 0.0
    return np.arctan2(x, y)

def segment_length(x1: float, y1: float, z1: float, 
                  x2: float, y2: float, z2: float) -> float:
    """Calculate segment length"""
    return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

def direction_cosines(x1: float, y1: float, z1: float,
                     x2: float, y2: float, z2: float) -> Tuple[float, float, float]:
    """Calculate direction cosines for a segment"""
    length = segment_length(x1, y1, z1, x2, y2, z2)
    if length == 0:
        return 0.0, 0.0, 0.0
    
    cab = (x2 - x1) / length  # Alpha
    sab = (y2 - y1) / length  # Beta
    salp = (z2 - z1) / length  # Gamma
    
    return cab, sab, salp

# ============================================================================
# PLACEHOLDER FUNCTIONS FOR UNIMPLEMENTED GEOMETRY TYPES
# ============================================================================

def process_wire_grid(nec_data: NECData, itg: int, ns: int,
                     x1: float, y1: float, z1: float,
                     x2: float, y2: float, z2: float, rad: float) -> None:
    """Process wire grid geometry (GR) - placeholder"""
    print(f"Wire Grid: Not yet implemented")

def process_wire_sphere(nec_data: NECData, itg: int, ns: int,
                       x1: float, y1: float, z1: float,
                       x2: float, y2: float, z2: float, rad: float) -> None:
    """Process wire sphere geometry (GS) - placeholder"""
    print(f"Wire Sphere: Not yet implemented")

def process_wire_ellipse(nec_data: NECData, itg: int, ns: int,
                        x1: float, y1: float, z1: float,
                        x2: float, y2: float, z2: float, rad: float) -> None:
    """Process wire ellipse geometry (GE) - placeholder"""
    print(f"Wire Ellipse: Not yet implemented")

def process_wire_move(nec_data: NECData, itg: int, ns: int,
                     x1: float, y1: float, z1: float,
                     x2: float, y2: float, z2: float, rad: float) -> None:
    """Process wire move geometry (GM) - placeholder"""
    print(f"Wire Move: Not yet implemented")

def process_surface_move(nec_data: NECData, itg: int, ns: int,
                        x1: float, y1: float, z1: float,
                        x2: float, y2: float, z2: float, rad: float) -> None:
    """Process surface move geometry (SM) - placeholder"""
    print(f"Surface Move: Not yet implemented")

def process_surface_file(nec_data: NECData) -> None:
    """Process surface file geometry (GF) - placeholder"""
    print(f"Surface File: Not yet implemented")

def process_surface_arc(nec_data: NECData, itg: int, ns: int,
                       x1: float, y1: float, z1: float,
                       x2: float, y2: float, z2: float, rad: float) -> None:
    """Process surface arc geometry (GA) - placeholder"""
    print(f"Surface Arc: Not yet implemented")

def process_surface_cylinder(nec_data: NECData, itg: int, ns: int,
                            x1: float, y1: float, z1: float,
                            x2: float, y2: float, z2: float, rad: float) -> None:
    """Process surface cylinder geometry (GC) - placeholder"""
    print(f"Surface Cylinder: Not yet implemented") 