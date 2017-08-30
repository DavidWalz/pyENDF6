"""
Tools to read an ENDF-6 data file.

From https://t2.lanl.gov/nis/endf/intro05.html
An ENDF-format nuclear data library has an hierarchical structure by tape, material, file, and section, denoted by numerical identifiers.

Tape is a data-file that contains one or more ENDF materials in increasing order by MAT.
    Each material contains several files in increasing order by MF.
    Each file contains several sections in increasing order by MT.

MAT labels an the target material as ZZXX, where XX starts from 25 for the lightest common isotope and increases in steps of 3 to allow for isomers.
    MAT=125   H-1
    MAT=128   H-2
    MAT=2625  Fe-54

MF labels an ENDF file to store different types of data:
    MF=1  descriptive and miscellaneous data,
    MF=2  resonance parameter data,
    MF=3  reaction cross sections vs energy,
    MF=4  angular distributions,
    MF=5  energy distributions,
    MF=6  energy-angle distributions,
    MF=7  thermal scattering data,
    MF=8  radioactivity data
    MF=9-10  nuclide production data,
    MF=12-15  photon production data, and
    MF=30-36  covariance data.

MT labels an ENDF section, usually used to hold different reactions, e.g.
    MT=1   total cross section
    MT=2   elastic scattering
    MT=3  Total photo-absorption cross section
    MT=16  (n,2n) reaction
    MT=18  fission
    MT=102 radiative capture
"""

import numpy as np


slices = {
    'MAT': slice(66, 70),
    'MF': slice(70, 72),
    'MT': slice(72, 75),
    'line': slice(75, 80),
    'content': slice(0, 66),
    'data': (slice(0, 11), slice(11, 22), slice(22, 33), slice(33, 44), slice(44, 55), slice(55, 66))}


def read_float(v):
    """
    Convert ENDF6 string to float
    (the ENDF6 float representation omits the e for exponent and may contain blanks)
    """
    try:
        number = float(v[0] + v[1:].replace(' ', ''))
    except ValueError:
        number = float(v[0] + v[1:].replace(' ', '').replace('+', 'e+').replace('-', 'e-'))
    return number


def read_line(l):
    """Read first 6*11 characters of a line as floats"""
    return [read_float(l[s]) for s in slices['data']]


def read_table(lines):
    """
    Parse tabulated data in a section
    https://t2.lanl.gov/nis/endf/intro07.html
    https://t2.lanl.gov/nis/endf/intro08.html
    https://t2.lanl.gov/nis/endf/intro09.html
    """
    # header line 1: (100*Z+A), mass in [m_neutron]
    # [MAT, 3, MT/ ZA, AWR, 0, 0, 0, 0] HEAD

    # header line 2: Q-value and some counts
    # [MAT, 3, MT/ QM, QI, 0, LR, NR, NP/ EINT/ S(E)] TAB1
    f = read_line(lines[1])
    nS = int(f[4])  # number of interpolation sections
    nP = int(f[5])  # number of data points

    # header line 3: interpolation information
    # [MAT, 3, 0/ 0.0, 0.0, 0, 0, 0, 0] SEND
    # 1   y is constant in x (constant, histogram)
    # 2   y is linear in x (linear-linear)
    # 3   y is linear in ln(x) (linear-log)
    # 4   ln(y) is linear in x (log-linear)
    # 5   ln(y) is linear in ln(x) (log-log)
    # 6   y obeys a Gamow charged-particle penetrability law

    # data lines
    x = []
    y = []
    for l in lines[3:-1]:
        f = read_line(l)
        x.append(f[0])
        y.append(f[1])
        x.append(f[2])
        y.append(f[3])
        x.append(f[4])
        y.append(f[5])
    return np.array(x[0:nP]), np.array(y[0:nP])


def find_file(lines, MF=1):
    """Locate and return a certain section"""
    v = [l[slices['MF']] for l in lines]
    n = len(v)
    cmpstr = '%2s' % MF       # search string
    i0 = v.index(cmpstr)            # first occurrence
    i1 = n - v[::-1].index(cmpstr)  # last occurrence
    return lines[i0: i1]


def find_section(lines, MF=3, MT=3):
    """Locate and return a certain section"""
    v = [l[70:75] for l in lines]
    n = len(v)
    cmpstr = '%2s%3s' % (MF, MT)       # search string
    i0 = v.index(cmpstr)            # first occurrence
    i1 = n - v[::-1].index(cmpstr)  # last occurrence
    return lines[i0: i1]


def list_content(lines):
    """Return set of unique tuples (MAT, MF, MT)"""
    s0 = slices['MAT']
    s1 = slices['MF']
    s2 = slices['MT']
    content = set(((int(l[s0]), int(l[s1]), int(l[s2])) for l in lines))

    # remove section delimiters
    for c in content.copy():
        if 0 in c:
            content.discard(c)
    return content
