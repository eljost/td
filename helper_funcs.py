#!/usr/bin/env python
# -*- coding: utf-8 -*-

def eVformat(eVs):
    fmt_str = "{:.2f}"
    return [fmt_str.format(eV) for eV in eVs]

def fformat(fs):
    fmt_str = "\\num{{{:.2e}}}"
    return [fmt_str.format(f) for f in fs]

def numformat(fs, places=2):
    fmt_str = "\\num{{{:.{places}f}}}"
    return [fmt_str.format(f, places=places) for f in fs]
