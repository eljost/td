#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path
import re

import numpy as np
from mytabulate import tabulate

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
