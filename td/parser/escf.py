#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import re

from td.helper_funcs import HARTREE2EV, HARTREE2NM
from td.ExcitedState import ExcitedState


def parse_escf(text):
    # In openshell calculations TURBOMOLE omits the multiplicity in
    # the string.
    sym = "(\d+)\s+(singlet|doublet|triplet|quartet|quintet|sextet)?" \
          "\s*([\w'\"]+)\s+excitation"
    sym_re = re.compile(sym)
    syms = sym_re.findall(text)
    syms = [(int(id_), spin, spat) for id_, spin, spat in syms]

    exc_energy = "Excitation energy:\s*([\d\.E\+-]+)"
    exc_energy_re = re.compile(exc_energy)
    ees = exc_energy_re.findall(text)
    ees = [float(ee) for ee in ees]

    osc_strength = "mixed representation:\s*([\d\.E\+-]+)"
    osc_strength_re = re.compile(osc_strength)
    oscs = osc_strength_re.findall(text)
    oscs = [float(osc) for osc in oscs]

    dom_contrib = "2\*100(.*?)Change of electron number"
    dom_contrib_re = re.compile(dom_contrib, flags=re.MULTILINE | re.DOTALL)
    dcs = dom_contrib_re.findall(text)
    dc_str = "(\d+) ([\w'\"]+)\s*(beta|alpha)?\s+([-\d\.]+)\s*" \
             "(\d+) ([\w'\"]+)\s*(beta|alpha)?\s+([-\d\.]+)\s*" \
             "([\d\.]+)"
    dc_re = re.compile(dc_str)
    dcs_parsed = [dc_re.findall(exc) for exc in dcs]

    excited_states = list()
    dom_contribs = list()
    for sym, ee, osc, dc in zip(syms, ees, oscs, dcs_parsed):
        id_, spin, spat = sym
        dE = ee * HARTREE2EV
        l = HARTREE2NM / ee

        exc_state = ExcitedState(id_, spin, spat, dE, l, osc, "???")
        excited_states.append(exc_state)
        for d in dc:
            start_mo = d[0]
            start_irrep = d[1]
            start_spin = d[2]
            final_mo = d[4]
            final_irrep = d[5]
            final_spin = d[6]
            to_or_from = "->"
            contrib = float(d[8]) / 100
            exc_state.add_mo_transition(start_mo, to_or_from, final_mo,
                                        ci_coeff=-0, contrib=contrib,
                                        start_spin=start_spin,
                                        final_spin=final_spin,
                                        start_irrep=start_irrep,
                                        final_irrep=final_irrep)

    return excited_states, dom_contribs
