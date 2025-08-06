"""
Main Program Module for Python NEC-2 Implementation

This module contains the main program logic converted from the Fortran code,
including the main loop, frequency sweeps, and program control flow.
"""

import numpy as np
from typing import Optional, List, Dict, Any
from .constants import *
from .data_structures import *
from .geometry import *
from .matrix_operations import *
from .field_calculations import *
from .ground_effects import *
from .network_analysis import *
from .input_output import *

# ============================================================================
# MAIN PROGRAM LOGIC
# ============================================================================

def main():
    """
    Main program entry point (converted from Fortran main program)
    """
    # Initialize main data structure
    nec_data = NECData()
    
    # Print program header
    print(format_header())
    
    # Main program loop
    while True:
        try:
            # Read control card
            card_type, params = read_control_card()
            
            if card_type == 'FR':
                # Frequency card
                process_frequency_card(params, nec_data)
            elif card_type == 'LD':
                # Loading card
                process_loading_card(params, nec_data)
            elif card_type == 'GN':
                # Ground card
                process_ground_card(params, nec_data)
            elif card_type == 'EX':
                # Excitation card
                process_excitation_card(params, nec_data)
            elif card_type == 'NT':
                # Network card
                process_network_card(params, nec_data)
            elif card_type == 'XQ':
                # Execute card
                execute_analysis(nec_data)
            elif card_type == 'NE':
                # Near field card
                process_near_field_card(params, nec_data)
            elif card_type == 'GD':
                # Ground representation card
                process_ground_representation_card(params, nec_data)
            elif card_type == 'RP':
                # Radiation pattern card
                process_radiation_pattern_card(params, nec_data)
            elif card_type == 'CP':
                # Coupling card
                process_coupling_card(params, nec_data)
            elif card_type == 'EN':
                # End card
                break
            elif card_type == 'CE' or card_type == 'CM':
                # Comment card
                print(f"COMMENT: {params}")
            else:
                print(f"Unknown card type: {card_type}")
                
        except EOFError:
            print("End of input file reached")
            break
        except Exception as e:
            print(f"Error processing card: {e}")
            break
    
    print("NEC-2 Analysis Complete")

def read_control_card() -> tuple:
    """
    Read a control card from input
    
    Returns:
        Tuple of (card_type, parameters)
    """
    # This is a simplified implementation
    # In the full implementation, this would read from the actual input file
    
    # Placeholder - would read actual card
    card_type = "FR"
    params = {"freq": 300.0, "nfrq": 1, "dfrq": 0.0, "freq2": 0.0}
    
    return card_type, params

def process_frequency_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process frequency card (FR from Fortran)
    
    Args:
        params: Frequency parameters
        nec_data: Main NEC data structure
    """
    freq = params.get('freq', DEFAULT_FREQ_MHZ)
    nfrq = params.get('nfrq', 1)
    dfrq = params.get('dfrq', 0.0)
    freq2 = params.get('freq2', 0.0)
    
    # Store frequency parameters
    nec_data.fmhz = freq
    nec_data.frequency.nfrq = nfrq
    nec_data.frequency.dfrq = dfrq
    nec_data.frequency.freq2 = freq2
    
    # Calculate wavelength
    nec_data.geometry.wlam = CVEL / freq
    
    print(f"FREQUENCY: {freq:.4e} MHZ")
    print(f"WAVELENGTH: {nec_data.geometry.wlam:.4e} METERS")

def process_loading_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process loading card (LD from Fortran)
    
    Args:
        params: Loading parameters
        nec_data: Main NEC data structure
    """
    ldtyp = params.get('ldtyp', 0)
    ldtag = params.get('ldtag', 0)
    ldtagf = params.get('ldtagf', 0)
    ldtagt = params.get('ldtagt', 0)
    zlr = params.get('zlr', 0.0)
    zli = params.get('zli', 0.0)
    zlc = params.get('zlc', 0.0)
    
    # Check if we have room for more loading cards
    if nec_data.loading.nload >= LOADMX:
        print("ERROR: NUMBER OF LOADING CARDS EXCEEDS STORAGE ALLOTTED")
        return
    
    # Store loading data
    idx = nec_data.loading.nload
    nec_data.loading.ldtyp[idx] = ldtyp
    nec_data.loading.ldtag[idx] = ldtag
    nec_data.loading.ldtagf[idx] = ldtagf
    nec_data.loading.ldtagt[idx] = ldtagt
    nec_data.loading.zlr[idx] = zlr
    nec_data.loading.zli[idx] = zli
    nec_data.loading.zlc[idx] = zlc
    
    nec_data.loading.nload += 1
    
    print(f"LOADING CARD {nec_data.loading.nload}: TYPE={ldtyp}, TAG={ldtag}")

def process_ground_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process ground card (GN from Fortran)
    
    Args:
        params: Ground parameters
        nec_data: Main NEC data structure
    """
    iperf = params.get('iperf', 0)
    nradl = params.get('nradl', 0)
    epsr = params.get('epsr', 1.0)
    sig = params.get('sig', 0.0)
    scrwlt = params.get('scrwlt', 0.0)
    scrwrt = params.get('scrwrt', 0.0)
    
    # Store ground parameters
    nec_data.ground.iperf = iperf
    nec_data.ground.nradl = nradl
    nec_data.epsr = epsr
    nec_data.sig = sig
    nec_data.scrwlt = scrwlt
    nec_data.scrwrt = scrwrt
    
    # Calculate ground parameters
    if iperf == 0:
        # Free space
        nec_data.ground.zrati = 1.0 + 0j
        nec_data.ground.zrati2 = 1.0 + 0j
        nec_data.ground.frati = 1.0 + 0j
        print("GROUND: FREE SPACE")
    elif iperf == 1:
        # Perfect ground
        nec_data.ground.zrati = 0.0 + 0j
        nec_data.ground.zrati2 = 0.0 + 0j
        nec_data.ground.frati = 1.0 + 0j
        print("GROUND: PERFECT GROUND")
    elif iperf == 2:
        # Finite ground (Sommerfeld)
        nec_data.ground.zrati = np.sqrt(1.0 / complex(epsr, -sig * nec_data.geometry.wlam * 59.96))
        nec_data.ground.zrati2 = nec_data.ground.zrati
        nec_data.ground.frati = 1.0 + 0j
        print(f"GROUND: FINITE GROUND (EPSR={epsr}, SIG={sig})")
    
    # Set symmetry flag
    if iperf > 0:
        nec_data.ground.ksymp = 2
    else:
        nec_data.ground.ksymp = 1

def process_excitation_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process excitation card (EX from Fortran)
    
    Args:
        params: Excitation parameters
        nec_data: Main NEC data structure
    """
    itype = params.get('itype', 0)
    itag = params.get('itag', 0)
    iseg = params.get('iseg', 0)
    i1 = params.get('i1', 0)
    i2 = params.get('i2', 0)
    f1 = params.get('f1', 0.0)
    f2 = params.get('f2', 0.0)
    f3 = params.get('f3', 0.0)
    f4 = params.get('f4', 0.0)
    f5 = params.get('f5', 0.0)
    f6 = params.get('f6', 0.0)
    
    # Check if we have room for more excitation cards
    if nec_data.excitation.nsant >= NSMAX:
        print("ERROR: NUMBER OF EXCITATION CARDS EXCEEDS STORAGE ALLOTTED")
        return
    
    # Store excitation data
    idx = nec_data.excitation.nsant
    nec_data.excitation.itype[idx] = itype
    nec_data.excitation.itag[idx] = itag
    nec_data.excitation.iseg[idx] = iseg
    nec_data.excitation.i1[idx] = i1
    nec_data.excitation.i2[idx] = i2
    nec_data.excitation.f1[idx] = f1
    nec_data.excitation.f2[idx] = f2
    nec_data.excitation.f3[idx] = f3
    nec_data.excitation.f4[idx] = f4
    nec_data.excitation.f5[idx] = f5
    nec_data.excitation.f6[idx] = f6
    
    # Calculate source parameters based on type
    if itype == 0:
        # Voltage source
        nec_data.excitation.vsant[idx] = complex(f1, f2)
        nec_data.excitation.isant[idx] = iseg
        print(f"EXCITATION: VOLTAGE SOURCE ON SEGMENT {iseg}")
    elif itype == 1:
        # Current source
        nec_data.excitation.vqd[idx] = complex(f1, f2)
        nec_data.excitation.ivqd[idx] = iseg
        print(f"EXCITATION: CURRENT SOURCE ON SEGMENT {iseg}")
    
    nec_data.excitation.nsant += 1

def process_network_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process network card (NT from Fortran)
    
    Args:
        params: Network parameters
        nec_data: Main NEC data structure
    """
    ntyp = params.get('ntyp', 1)
    iseg1 = params.get('iseg1', 0)
    iseg2 = params.get('iseg2', 0)
    x11r = params.get('x11r', 0.0)
    x11i = params.get('x11i', 0.0)
    x12r = params.get('x12r', 0.0)
    x12i = params.get('x12i', 0.0)
    x22r = params.get('x22r', 0.0)
    x22i = params.get('x22i', 0.0)
    
    # Check if we have room for more network cards
    if nec_data.network.nonet >= NETMX:
        print("ERROR: NUMBER OF NETWORK CARDS EXCEEDS STORAGE ALLOTTED")
        return
    
    # Store network data
    idx = nec_data.network.nonet
    nec_data.network.ntyp[idx] = ntyp
    nec_data.network.iseg1[idx] = iseg1
    nec_data.network.iseg2[idx] = iseg2
    nec_data.network.x11r[idx] = x11r
    nec_data.network.x11i[idx] = x11i
    nec_data.network.x12r[idx] = x12r
    nec_data.network.x12i[idx] = x12i
    nec_data.network.x22r[idx] = x22r
    nec_data.network.x22i[idx] = x22i
    
    nec_data.network.nonet += 1
    
    network_type = NETWORK_TYPES.get(ntyp, "UNKNOWN")
    print(f"NETWORK: {network_type} BETWEEN SEGMENTS {iseg1} AND {iseg2}")

def process_near_field_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process near field card (NE from Fortran)
    
    Args:
        params: Near field parameters
        nec_data: Main NEC data structure
    """
    xnr = params.get('xnr', 0.0)
    ynr = params.get('ynr', 0.0)
    znr = params.get('znr', 0.0)
    dxnr = params.get('dxnr', 0.0)
    dynr = params.get('dynr', 0.0)
    dznr = params.get('dznr', 0.0)
    nrx = params.get('nrx', 1)
    nry = params.get('nry', 1)
    nrz = params.get('nrz', 1)
    
    # Store near field parameters
    nec_data.pattern.xnr = xnr
    nec_data.pattern.ynr = ynr
    nec_data.pattern.znr = znr
    nec_data.pattern.dxnr = dxnr
    nec_data.pattern.dynr = dynr
    nec_data.pattern.dznr = dznr
    nec_data.pattern.nrx = nrx
    nec_data.pattern.nry = nry
    nec_data.pattern.nrz = nrz
    nec_data.pattern.near = 1
    
    print(f"NEAR FIELD: {nrx}x{nry}x{nrz} POINTS STARTING AT ({xnr:.3f}, {ynr:.3f}, {znr:.3f})")

def process_ground_representation_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process ground representation card (GD from Fortran)
    
    Args:
        params: Ground representation parameters
        nec_data: Main NEC data structure
    """
    ifar = params.get('ifar', 0)
    nradl = params.get('nradl', 0)
    epsr2 = params.get('epsr2', 1.0)
    sig2 = params.get('sig2', 0.0)
    clt = params.get('clt', 0.0)
    cht = params.get('cht', 0.0)
    
    # Store ground representation parameters
    nec_data.ground.ifar = ifar
    nec_data.ground.nradl = nradl
    nec_data.pattern.epsr2 = epsr2
    nec_data.pattern.sig2 = sig2
    nec_data.pattern.clt = clt
    nec_data.pattern.cht = cht
    
    print(f"GROUND REPRESENTATION: IFAR={ifar}, NRADL={nradl}")

def process_radiation_pattern_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process radiation pattern card (RP from Fortran)
    
    Args:
        params: Radiation pattern parameters
        nec_data: Main NEC data structure
    """
    ipd = params.get('ipd', 0)
    iavp = params.get('iavp', 0)
    inor = params.get('inor', 0)
    iax = params.get('iax', 0)
    ixtyp = params.get('ixtyp', 0)
    nth = params.get('nth', 91)
    nph = params.get('nph', 1)
    thets = params.get('thets', 0.0)
    phis = params.get('phis', 0.0)
    dth = params.get('dth', 1.0)
    dph = params.get('dph', 0.0)
    rfld = params.get('rfld', 0.0)
    gnor = params.get('gnor', 1.0)
    
    # Store radiation pattern parameters
    nec_data.pattern.ipd = ipd
    nec_data.pattern.iavp = iavp
    nec_data.pattern.inor = inor
    nec_data.pattern.iax = iax
    nec_data.pattern.ixtyp = ixtyp
    nec_data.pattern.nth = nth
    nec_data.pattern.nph = nph
    nec_data.pattern.thets = thets
    nec_data.pattern.phis = phis
    nec_data.pattern.dth = dth
    nec_data.pattern.dph = dph
    nec_data.pattern.rfld = rfld
    nec_data.pattern.gnor = gnor
    nec_data.pattern.ifar = 1
    
    print(f"RADIATION PATTERN: {nth}x{nph} POINTS FROM ({thets:.1f}°, {phis:.1f}°)")

def process_coupling_card(params: Dict[str, Any], nec_data: NECData) -> None:
    """
    Process coupling card (CP from Fortran)
    
    Args:
        params: Coupling parameters
        nec_data: Main NEC data structure
    """
    ncoup = params.get('ncoup', 0)
    nctag = params.get('nctag', [0] * 5)
    ncseg = params.get('ncseg', [0] * 5)
    
    # Check if we have room for coupling calculation
    if ncoup > 5:
        print("ERROR: NUMBER OF SEGMENTS IN COUPLING CALCULATION (CP) EXCEEDS LIMIT")
        return
    
    # Store coupling parameters
    nec_data.coupling.ncoup = ncoup
    for i in range(ncoup):
        nec_data.coupling.nctag[i] = nctag[i]
        nec_data.coupling.ncseg[i] = ncseg[i]
    
    print(f"COUPLING: {ncoup} SEGMENTS")

def execute_analysis(nec_data: NECData) -> None:
    """
    Execute the main analysis (XQ from Fortran)
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\n" + "="*60)
    print("EXECUTING NEC-2 ANALYSIS")
    print("="*60)
    
    # Frequency loop
    for ifrq in range(nec_data.frequency.nfrq):
        # Calculate frequency for this iteration
        if nec_data.frequency.nfrq == 1:
            freq = nec_data.fmhz
        else:
            freq = nec_data.fmhz + ifrq * nec_data.frequency.dfrq
        
        # Update wavelength
        nec_data.geometry.wlam = CVEL / freq
        
        print(f"\nFREQUENCY {ifrq + 1}: {freq:.4e} MHZ")
        print(f"WAVELENGTH: {nec_data.geometry.wlam:.4e} METERS")
        
        # Generate geometry data
        datagn(nec_data)
        
        # Apply loading
        apply_loading(nec_data)
        
        # Set up ground effects
        setup_ground_effects(nec_data)
        
        # Set up interaction matrix
        setup_matrix(nec_data)
        
        # Factor matrix
        factor_matrix(nec_data)
        
        # Solve for currents
        solve_currents(nec_data)
        
        # Calculate fields
        calculate_fields(nec_data)
        
        # Print results
        print_results(nec_data, ifrq)

def datagn(nec_data: NECData) -> None:
    """
    Generate geometry data (DATAGN from Fortran)
    
    Args:
        nec_data: Main NEC data structure
    """
    print("\nGENERATING GEOMETRY DATA...")
    
    # Read geometry cards until GE card
    while True:
        card_type, params = read_geometry_card()
        
        if card_type == 'GE':
            # End of geometry
            break
        elif card_type == 'GW':
            # Wire
            process_wire(nec_data, params['nwire'], params['itg'], params['ns'],
                        params['xw1'], params['yw1'], params['zw1'],
                        params['xw2'], params['yw2'], params['zw2'], params['rad'])
        elif card_type == 'GA':
            # Wire arc
            process_wire_arc(nec_data, params['nwire'], params['itg'], params['ns'],
                           params['xc'], params['yc'], params['zc'], params['radius'])
        elif card_type == 'GH':
            # Helix
            process_helix(nec_data, params['nwire'], params['itg'], params['ns'],
                         params['xc'], params['yc'], params['zc'],
                         params['xh'], params['yh'], params['zh'], params['rad'])
        elif card_type == 'SP':
            # Surface patch
            process_surface_patch(nec_data, params['itg'], params['ns'],
                                params['x1'], params['y1'], params['z1'],
                                params['x2'], params['y2'], params['z2'], params['rad'])
        else:
            print(f"Unknown geometry card: {card_type}")
    
    # Set up connections
    setup_connections(nec_data)
    
    print(f"GEOMETRY COMPLETE: {nec_data.geometry.n} SEGMENTS, {nec_data.geometry.m} PATCHES")

def apply_loading(nec_data: NECData) -> None:
    """
    Apply loading to segments
    
    Args:
        nec_data: Main NEC data structure
    """
    if nec_data.loading.nload == 0:
        return
    
    print("APPLYING LOADING...")
    
    for i in range(nec_data.loading.nload):
        ldtyp = nec_data.loading.ldtyp[i]
        ldtag = nec_data.loading.ldtag[i]
        ldtagf = nec_data.loading.ldtagf[i]
        ldtagt = nec_data.loading.ldtagt[i]
        zlr = nec_data.loading.zlr[i]
        zli = nec_data.loading.zli[i]
        zlc = nec_data.loading.zlc[i]
        
        # Apply loading to segments
        for j in range(nec_data.geometry.n):
            if nec_data.geometry.itag[j] == ldtag:
                if ldtagf == 0 or ldtagf <= j + 1 <= ldtagt:
                    zload = complex(zlr, zli) * zlc
                    nec_data.loading.zarray[j] = zload

def setup_ground_effects(nec_data: NECData) -> None:
    """
    Set up ground effects
    
    Args:
        nec_data: Main NEC data structure
    """
    if nec_data.ground.iperf == 2:
        # Sommerfeld ground
        print("SETTING UP SOMMERFELD GROUND...")
        som2d(nec_data.fmhz, nec_data.epsr, nec_data.sig)

def setup_matrix(nec_data: NECData) -> None:
    """
    Set up interaction matrix
    
    Args:
        nec_data: Main NEC data structure
    """
    print("SETTING UP INTERACTION MATRIX...")
    
    neq = nec_data.geometry.n + 2 * nec_data.geometry.m
    npeq = nec_data.geometry.np + 2 * nec_data.geometry.mp
    
    # Initialize matrix
    cm = np.zeros((neq, npeq), dtype=np.complex128)
    
    # Fill matrix
    cmset(nec_data, neq, cm, DEFAULT_RKH, 0)
    
    # Store matrix
    nec_data.matrix.cm = cm
    nec_data.matrix.icase = 1
    nec_data.matrix.imat = neq * npeq

def factor_matrix(nec_data: NECData) -> None:
    """
    Factor the interaction matrix
    
    Args:
        nec_data: Main NEC data structure
    """
    print("FACTORING MATRIX...")
    
    neq = nec_data.geometry.n + 2 * nec_data.geometry.m
    npeq = nec_data.geometry.np + 2 * nec_data.geometry.mp
    
    # Initialize pivot array
    ip = np.zeros(neq, dtype=np.int32)
    
    # Factor matrix
    factrs(npeq, neq, nec_data.matrix.cm, ip, nec_data.ip, 11, 12, 13, 14)
    
    # Store pivot array
    nec_data.ip = ip

def solve_currents(nec_data: NECData) -> None:
    """
    Solve for currents
    
    Args:
        nec_data: Main NEC data structure
    """
    print("SOLVING FOR CURRENTS...")
    
    neq = nec_data.geometry.n + 2 * nec_data.geometry.m
    npeq = nec_data.geometry.np + 2 * nec_data.geometry.mp
    
    # Initialize right-hand side
    b = np.zeros(neq, dtype=np.complex128)
    
    # Set up excitation
    for i in range(nec_data.excitation.nsant):
        iseg = nec_data.excitation.isant[i] - 1
        vs = nec_data.excitation.vsant[i]
        b[iseg] = vs
    
    # Solve system
    solve(npeq, neq, nec_data.matrix.cm, nec_data.ip, nec_data.ip, 11, 12, 13, 14)
    
    # Store currents
    nec_data.current.cur[:neq] = b[:neq]

def calculate_fields(nec_data: NECData) -> None:
    """
    Calculate fields
    
    Args:
        nec_data: Main NEC data structure
    """
    print("CALCULATING FIELDS...")
    
    # Calculate near field if requested
    if nec_data.pattern.near > 0:
        nfpat(nec_data)
    
    # Calculate radiation pattern if requested
    if nec_data.pattern.ifar > 0:
        rdpat(nec_data)
    
    # Calculate coupling if requested
    if nec_data.coupling.ncoup > 0:
        couple(nec_data.current.cur, nec_data.geometry.wlam, nec_data)

def print_results(nec_data: NECData, ifrq: int) -> None:
    """
    Print analysis results
    
    Args:
        nec_data: Main NEC data structure
        ifrq: Frequency index
    """
    print(f"\nRESULTS FOR FREQUENCY {ifrq + 1}:")
    print("-" * 40)
    
    # Print current summary
    print("CURRENT SUMMARY:")
    for i in range(min(10, nec_data.geometry.n)):
        cur_mag = abs(nec_data.current.cur[i])
        cur_phase = np.angle(nec_data.current.cur[i], deg=True)
        print(f"  SEG {i+1}: {cur_mag:.3e} A, {cur_phase:.1f}°")
    
    # Print input impedance if available
    if nec_data.excitation.nsant > 0:
        iseg = nec_data.excitation.isant[0] - 1
        vs = nec_data.excitation.vsant[0]
        cur = nec_data.current.cur[iseg]
        if abs(cur) > 1e-20:
            zin = vs / cur
            print(f"INPUT IMPEDANCE: {zin:.3e} OHMS")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def read_geometry_card() -> tuple:
    """
    Read a geometry card from input
    
    Returns:
        Tuple of (card_type, parameters)
    """
    # This is a simplified implementation
    # In the full implementation, this would read from the actual input file
    
    # Placeholder - would read actual geometry card
    card_type = "GE"
    params = {}
    
    return card_type, params

if __name__ == "__main__":
    main() 