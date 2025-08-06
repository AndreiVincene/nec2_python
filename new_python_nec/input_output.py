"""
Input/Output Module for Python NEC-2 Implementation

This module contains all input/output functions converted from the Fortran code,
including READMN, READGM, PARSIT, and other I/O subroutines.
"""

import numpy as np
from typing import Tuple, List, Optional, Dict, Any
from .constants import *
from .data_structures import *

# ============================================================================
# MAIN INPUT ROUTINE (READMN)
# ============================================================================

def readmn(inunit: int, code: str, i1: int, i2: int, i3: int, i4: int,
           f1: float, f2: float, f3: float, f4: float, f5: float, f6: float) -> None:
    """
    Read a control record and parse it (READMN from Fortran)
    
    Args:
        inunit: Input unit number
        code: Two letter mnemonic code
        i1-i4: Integer values from record
        f1-f6: Real values from record
    """
    intval = [0] * 4
    reaval = [0.0] * 6
    
    # Call the routine to read the record and parse it
    iEOF = parsit(inunit, 4, 6, code, intval, reaval)
    
    # Set the return variables to the buffer array elements
    if iEOF < 0:
        code = 'EN'
    
    i1 = intval[0]
    i2 = intval[1]
    i3 = intval[2]
    i4 = intval[3]
    f1 = reaval[0]
    f2 = reaval[1]
    f3 = reaval[2]
    f4 = reaval[3]
    f5 = reaval[4]
    f6 = reaval[5]

# ============================================================================
# GEOMETRY INPUT ROUTINE (READGM)
# ============================================================================

def readgm(inunit: int, code: str, i1: int, i2: int,
           r1: float, r2: float, r3: float, r4: float, 
           r5: float, r6: float, r7: float) -> None:
    """
    Read a geometry record and parse it (READGM from Fortran)
    
    Args:
        inunit: Input unit number
        code: Two letter mnemonic code
        i1, i2: Integer values from record
        r1-r7: Real values from record
    """
    intval = [0] * 2
    reaval = [0.0] * 7
    
    # Call the routine to read the record and parse it
    iEOF = parsit(inunit, 2, 7, code, intval, reaval)
    
    # Set the return variables to the buffer array elements
    if iEOF < 0:
        code = 'GE'
    
    i1 = intval[0]
    i2 = intval[1]
    r1 = reaval[0]
    r2 = reaval[1]
    r3 = reaval[2]
    r4 = reaval[3]
    r5 = reaval[4]
    r6 = reaval[5]
    r7 = reaval[6]

# ============================================================================
# RECORD PARSING (PARSIT)
# ============================================================================

def parsit(inunit: int, maxint: int, maxrea: int, cmnd: str, 
           intfld: List[int], reafld: List[float]) -> int:
    """
    Read an input record and parse it (PARSIT from Fortran)
    
    Args:
        inunit: Input unit number
        maxint: Total number of integers in record
        maxrea: Total number of real values in record
        cmnd: Two letter mnemonic code
        intfld: Integer values from record
        reafld: Real values from record
        
    Returns:
        EOF status (negative if EOF)
    """
    # Global variables
    global ngfnam
    ngfnam = 'NGF2D.NEC'  # Default NGF filename
    
    # Read record
    try:
        rec = read_line(inunit)
        if rec is None:
            return -1
    except EOFError:
        return -1
    
    # Convert to upper case
    rec, totcol = upcase(rec)
    
    # Store opcode and clear field arrays
    cmnd = rec[0:2]
    intfld[:] = [0] * maxint
    reafld[:] = [0.0] * maxrea
    
    # Find field boundaries
    bgnfld = [0] * 12
    endfld = [0] * 12
    totfld = 0
    fldtrm = False
    
    # Parse fields
    for j in range(2, totcol):
        k = ord(rec[j])
        
        # Check for end of line comment ('!')
        if k == 33:  # '!'
            if fldtrm:
                endfld[totfld] = j - 1
            break
        
        # Set ending index for comma or space
        elif k == 32 or k == 44:  # space or comma
            if fldtrm:
                endfld[totfld] = j - 1
                fldtrm = False
        
        # Set beginning index for non-delimiter
        elif not fldtrm:
            totfld += 1
            fldtrm = True
            bgnfld[totfld] = j
    
    if fldtrm:
        endfld[totfld] = totcol
    
    # Handle special cases for WG and GF cards
    if cmnd in ['WG', 'GF']:
        ngfnam = 'NGF2D.NEC'  # Initialize default NGFNAM
    
    # Check field count
    if totfld == 0:
        return 0
    elif totfld > maxint + maxrea:
        print("***** CARD ERROR - TOO MANY FIELDS IN RECORD")
        print(f"***** TEXT -->  {rec}")
        raise ValueError("Card parsing error")
    
    # Parse integer values
    j = min(totfld, maxint)
    for i in range(j):
        length = endfld[i] - bgnfld[i] + 1
        buffer = rec[bgnfld[i]:endfld[i]]
        
        # Handle text field for WG/GF cards
        if cmnd in ['WG', 'GF'] and buffer[0] not in ['0', '1']:
            ngfnam = rec[bgnfld[i]:endfld[i]]
            return 0
        
        # Check for decimal point
        ind = buffer.find('.')
        if ind > 0 and ind < length:
            print(f"***** CARD ERROR - INVALID NUMBER AT INTEGER POSITION {i+1}")
            print(f"***** TEXT -->  {rec}")
            raise ValueError("Card parsing error")
        if ind == length:
            length -= 1
        
        try:
            intfld[i] = int(buffer[:length])
        except ValueError:
            print(f"***** CARD ERROR - INVALID NUMBER AT INTEGER POSITION {i+1}")
            print(f"***** TEXT -->  {rec}")
            raise ValueError("Card parsing error")
    
    # Parse real values
    if totfld > maxint:
        j = maxint + 1
        for i in range(j, totfld + 1):
            length = endfld[i] - bgnfld[i] + 1
            buffer = rec[bgnfld[i]:endfld[i]]
            
            # Add decimal point if missing
            ind = buffer.find('.')
            if ind == -1:
                inde = buffer.find('E')
                length += 1
                if inde == -1:
                    buffer += '.'
                else:
                    buffer = buffer[:inde] + '.' + buffer[inde:length-1]
            
            try:
                reafld[i - maxint] = float(buffer[:length])
            except ValueError:
                print(f"***** CARD ERROR - INVALID NUMBER AT REAL POSITION {i - maxint}")
                print(f"***** TEXT -->  {rec}")
                raise ValueError("Card parsing error")
    
    return 0

# ============================================================================
# CASE CONVERSION (UPCASE)
# ============================================================================

def upcase(intext: str) -> Tuple[str, int]:
    """
    Find the length of INTEXT and convert it to upper case (UPCASE from Fortran)
    
    Args:
        intext: Input text
        
    Returns:
        Tuple of (uppercase text, length)
    """
    length = len(intext)
    outtxt = ""
    
    for i in range(length):
        j = ord(intext[i])
        if j >= 96:
            j -= 32
        outtxt += chr(j)
    
    return outtxt, length

# ============================================================================
# FILE OPERATIONS
# ============================================================================

def read_line(inunit: int) -> Optional[str]:
    """
    Read a line from input unit
    
    Args:
        inunit: Input unit number
        
    Returns:
        Line content or None if EOF
    """
    # This is a simplified implementation
    # In the full implementation, this would read from the actual file
    try:
        # Placeholder for file reading
        return None
    except EOFError:
        return None

def write_line(outunit: int, line: str) -> None:
    """
    Write a line to output unit
    
    Args:
        outunit: Output unit number
        line: Line to write
    """
    # This is a simplified implementation
    # In the full implementation, this would write to the actual file
    print(line)

# ============================================================================
# PRINTING ROUTINES
# ============================================================================

def prnt(in1: int, in2: int, in3: int, fl1: float, fl2: float, fl3: float,
         fl4: float, fl5: float, fl6: float, ctype: str) -> None:
    """
    Print input data for impedance loading (PRNT from Fortran)
    
    Args:
        in1-in3: Integer values to be printed
        fl1-fl6: Real values to be printed
        ctype: Character string to be printed
    """
    cint = ['     '] * 3
    cflt = ['     '] * 6
    
    # Handle integer values
    if in1 == 0 and in2 == 0 and in3 == 0:
        cint[0] = '  ALL'
    else:
        if in1 != 0:
            cint[0] = f"{in1:5d}"
        if in2 != 0:
            cint[1] = f"{in2:5d}"
        if in3 != 0:
            cint[2] = f"{in3:5d}"
    
    # Handle real values
    if abs(fl1) > 1e-30:
        cflt[0] = f"{fl1:13.4e}"
    if abs(fl2) > 1e-30:
        cflt[1] = f"{fl2:13.4e}"
    if abs(fl3) > 1e-30:
        cflt[2] = f"{fl3:13.4e}"
    if abs(fl4) > 1e-30:
        cflt[3] = f"{fl4:13.4e}"
    if abs(fl5) > 1e-30:
        cflt[4] = f"{fl5:13.4e}"
    if abs(fl6) > 1e-30:
        cflt[5] = f"{fl6:13.4e}"
    
    # Print formatted line
    line = f"\n   {''.join(cint)}   {''.join(cflt)}   {ctype}"
    print(line)

# ============================================================================
# FORMAT STRINGS
# ============================================================================

def format_header() -> str:
    """Format program header"""
    return """*********************************************
NUMERICAL ELECTROMAGNETICS CODE (NEC-2D)
*********************************************"""

def format_frequency(freq: float, wlam: float) -> str:
    """Format frequency information"""
    return f"FREQUENCY={freq:.4e} MHZ\nWAVELENGTH={wlam:.4e} METERS"

def format_matrix_timing(fill_time: float, factor_time: float) -> str:
    """Format matrix timing information"""
    return f"FILL={fill_time:.3f} SEC., FACTOR={factor_time:.3f} SEC."

def format_power_budget(pin: float, prad: float, ploss: float, pnls: float, efficiency: float) -> str:
    """Format power budget information"""
    return f"""INPUT POWER={pin:.4e} WATTS
RADIATED POWER={prad:.4e} WATTS
STRUCTURE LOSS={ploss:.4e} WATTS
NETWORK LOSS={pnls:.4e} WATTS
EFFICIENCY={efficiency:.2f} PERCENT"""

def format_current_header() -> str:
    """Format current header"""
    return """SEG. TAG COORD. OF SEG. CENTER SEG. - - - CURRENT (AMPS) - - -
NO. NO. X Y Z LENGTH REAL IMAG. MAG. PHASE"""

def format_impedance_header() -> str:
    """Format impedance header"""
    return """FREQ. - UNNORMALIZED IMPEDANCE - - NORMALIZED IMPEDANCE -
MHZ RESISTANCE REACTANCE MAGNITUDE PHASE RESISTANCE REACTANCE MAGNITUDE PHASE"""

# ============================================================================
# ERROR HANDLING
# ============================================================================

def print_error(message: str, card_text: str = "") -> None:
    """
    Print error message
    
    Args:
        message: Error message
        card_text: Card text that caused the error
    """
    print(f"\n***** {message}")
    if card_text:
        print(f"***** TEXT -->  {card_text}")

def handle_card_error(error_type: str, position: int = 0, card_text: str = "") -> None:
    """
    Handle card parsing errors
    
    Args:
        error_type: Type of error
        position: Position where error occurred
        card_text: Card text that caused the error
    """
    if error_type == "TOO_MANY_FIELDS":
        print_error("CARD ERROR - TOO MANY FIELDS IN RECORD", card_text)
    elif error_type == "INVALID_INTEGER":
        print_error(f"CARD ERROR - INVALID NUMBER AT INTEGER POSITION {position}", card_text)
    elif error_type == "INVALID_REAL":
        print_error(f"CARD ERROR - INVALID NUMBER AT REAL POSITION {position}", card_text)
    else:
        print_error("CARD ERROR", card_text)
    
    raise ValueError("Card parsing error")

# ============================================================================
# FILE UNIT MANAGEMENT
# ============================================================================

class FileManager:
    """Manage file units and operations"""
    
    def __init__(self):
        self.files = {}
        self.ngfnam = 'NGF2D.NEC'
    
    def open_file(self, unit: int, filename: str, mode: str = 'r') -> None:
        """Open a file for the specified unit"""
        try:
            if mode == 'r':
                file_obj = open(filename, 'r')
            elif mode == 'w':
                file_obj = open(filename, 'w')
            elif mode == 'a':
                file_obj = open(filename, 'a')
            else:
                file_obj = open(filename, mode)
            
            self.files[unit] = file_obj
        except Exception as e:
            print(f"Error opening file {filename} for unit {unit}: {e}")
            raise
    
    def close_file(self, unit: int) -> None:
        """Close a file for the specified unit"""
        if unit in self.files:
            self.files[unit].close()
            del self.files[unit]
    
    def read_line(self, unit: int) -> Optional[str]:
        """Read a line from the specified unit"""
        if unit in self.files:
            try:
                return self.files[unit].readline().strip()
            except EOFError:
                return None
        return None
    
    def write_line(self, unit: int, line: str) -> None:
        """Write a line to the specified unit"""
        if unit in self.files:
            self.files[unit].write(line + '\n')
        else:
            print(line)  # Default to console output
    
    def set_ngf_filename(self, filename: str) -> None:
        """Set the NGF filename"""
        self.ngfnam = filename
    
    def get_ngf_filename(self) -> str:
        """Get the NGF filename"""
        return self.ngfnam

# Global file manager instance
file_manager = FileManager()

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_card_line(line: str) -> Tuple[str, List[int], List[float]]:
    """
    Parse a card line into components
    
    Args:
        line: Input line
        
    Returns:
        Tuple of (command, integers, reals)
    """
    # Remove comments
    if '!' in line:
        line = line[:line.index('!')]
    
    # Split into fields
    fields = line.split()
    if len(fields) < 2:
        return '', [], []
    
    command = fields[0]
    values = fields[1:]
    
    integers = []
    reals = []
    
    for value in values:
        try:
            # Try to parse as integer first
            int_val = int(value)
            integers.append(int_val)
        except ValueError:
            try:
                # Try to parse as float
                float_val = float(value)
                reals.append(float_val)
            except ValueError:
                # Skip invalid values
                continue
    
    return command, integers, reals

def format_output_line(values: List[Any], format_spec: str) -> str:
    """
    Format output line with specified format
    
    Args:
        values: Values to format
        format_spec: Format specification
        
    Returns:
        Formatted line
    """
    try:
        return format_spec.format(*values)
    except (ValueError, TypeError):
        # Fallback to simple formatting
        return ' '.join(str(v) for v in values) 