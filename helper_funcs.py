#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path
import re

import numpy as np
from mytabulate import tabulate

def nm2eV_ticks(wavelengths):
    return ["{:.2f}".format(nm2eV(nm)) for nm in wavelengths]

def get_items_with_biggest_gap(lst):
    gaps = [lst[i+1] - lst[i] for i in range(len(lst)-1)]
    index = np.argmax(gaps)
    return lst[index], lst[index+1]

def find_biggest_gap(lst):
    #if len(lst) is 1:
    #    return 0, lst[0], lst[0]
    gaps = [lst[i+1] - lst[i] for i in range(len(lst)-1)]
    index = np.argmax(gaps)
    return index, gaps[index], lst[index], lst[index+1]

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def swap_chars_in_str(strr, a, b):
    as_list = list(strr)
    as_list[a], as_list[b] = as_list[b], as_list[a]
    return "".join(as_list)

def hartree2eV(hartree):
    return hartree * 27.2114

def hartree2nm(hartree):
    return 6.626-34 * 2.9979e8 * 1e9 / (hartree * 27.2114 * 1.6e-19)

def nm2eV(nm):
    return 6.626e-34 / 1.6021e-19 * 2.9979e8 / (nm / 1.e9)

def eVformat(eVs):
    fmt_str = "{:.2f}"
    return [fmt_str.format(eV) for eV in eVs]

def fformat(fs):
    fmt_str = "\\num{{{:.2e}}}"
    return [fmt_str.format(f) for f in fs]

def numformat(fs, places=2):
    fmt_str = "\\num{{{:.{places}f}}}"
    return [fmt_str.format(f, places=places) for f in fs]

def to_spherical(cartesian):
    x, y, z = cartesian
    r = np.linalg.norm(cartesian)
    if np.isclose(r, 0):
        return (0, 0, 0)
    theta = np.arccos(z/r)
    if np.isclose(theta, 0):
        return (0, 0, r)
    elif np.isclose(theta, np.pi):
        return (0, 0, -r)
    phi = np.arctan2(y, x)
    return (r, theta, phi)

def print_deg(spherical):
    r, t, p = spherical
    print(r, np.rad2deg(t), np.rad2deg(p))

def count_core_electrons(fchk):
    atoms = fchk["Atomic numbers"]
    # A negative charges "adds" electrons...
    electrons = sum(atoms) - fchk["Charge"]
    # Calculate how many electrons are described by an ECP
    ecp_electrons = electrons - (fchk["Number of alpha electrons"] +
        fchk["Number of beta electrons"])

    core_electrons = 0
    for atom in atoms:
        if 3 <= atom <= 10:
            core_electrons += 2
            continue
        if 11 <= atom <= 36:
            core_electrons += 10
            continue
        if atom > 36:
            raise Exception("Calculation of core elctrons for this" \
                " atomic number not yet implemented.")
    core_electrons -= ecp_electrons

    return core_electrons

def split_in_occ_virt(chains):
    length = len(chains)
    assert((length % 2) == 0)

    occ_chains = chains[:length/2] 
    virt_chains = chains[length/2:]

    return occ_chains, virt_chains

def get_smallest_gap(chains):
    """Find smallest gap between chains representing occupied
    MOs and chains representing virtual MOs, usually the gap
    between the highest HOMO and the lowest LUMO."""

    occ_chains, virt_chains = split_in_occ_virt(chains)

    return np.min(virt_chains) - np.max(occ_chains)

def scan_variable(start, increment, steps):
    # Step 1 equals the starting geometry so always increment
    # the start value by  (i-1) times.
    return [start + increment * (i - 1) for i in steps]

def steps_from_fns(fns):
    return [int(re.search("step_(\d+)", fn).groups()[0]) for fn in fns]

def check_if_exists_else_get_handle(path, mode):
    if not os.path.exists(path):
        return open(path, mode)
    else:
        raise Exception("File '{0}' already exists".format(path))

def print_bt_table(output):
    # #, Sym., E in eV, E in nm, f, Typ, NTO, %
    # Expected argument format
    # Iterable of iterables with format
    # (Sym., E in eV, E in nm, f)

    # Replace E in eV and E in nm
    formatted_output = list()
    symms, EeVs, Enms, fs = zip(*output)
    EeVs = numformat(EeVs)
    Enms = numformat(Enms, 1)
    fs = fformat(fs)
    formatted_output = zip(
        range(1, len(symms)+1), symms, EeVs, Enms, fs)

    print tabulate(formatted_output, tablefmt="latex_booktabs")

if __name__ == "__main__":
    """
    print to_spherical((0,0,0))
    d = lambda s: print_deg(to_spherical(s))
    d((1,1,1))
    d((1,0,0))
    d((0,1,0))
    print to_spherical((0,0,1))
    d((0,0,1))
    """
    g = (1,2,3,4,5,8,13,29,16,15,14,88)
    print g
    print find_biggest_gap(g)
