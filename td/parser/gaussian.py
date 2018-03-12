#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

from td.ExcitedState import ExcitedState
from td.helper_funcs import conv

# Excitation energies and oscillator strengths
# Groups: id, spin symm, spat symm, dE, wavelength, osc. strength, S**2
EXC_LINE = "Excited State\s+(\d+):\s+([\w\.]+)-([\?\w'\"]+)\s+([0-9\.-]+) eV" \
           "\s+([0-9\.-]+) nm\s+f=([0-9\.-]+)\s+<S\*\*2>=([0-9\.]+)"
# ExcitedStates between MOs and corresponding CI-coefficients
# Groups: initial MO, final MO, ci coeffcient
TRS_LINE = r"([\dAB]+)\s*(->|<-)\s*([\dAB]+)\s*\s+([0-9\.-]+)"


def parse_tddft(text):
    # Determine multiplicity
    charge_mult_re = "\s*".join("Charge = ([\+\-\d]+) Multiplicity = (\d+)".split())
    _, mult = re.search(charge_mult_re, text).groups()
    mult = int(mult)
    assert(mult >= 1)

    lines = text.split("\n")
    excited_states = list()

    matched_exc_state = False
    for line in lines:
        line = line.strip()
        # look for initial and final MO and corresponding CI-coefficient
        if matched_exc_state:
            match_obj = re.match(TRS_LINE, line)
            if match_obj:
                group_list = list(match_obj.groups())
                start_mo, to_or_from, final_mo, ci_coeff = conv(group_list,
                                                                "sssf")
                def handle_mo_item(mo_item):
                    _, mo_num, spin = re.split("(\d+)", mo_item)
                    if not spin:
                        spin = "a"
                    return mo_num, spin
                start_mo = start_mo.lower()
                final_mo = final_mo.lower()
                start_mo_num, start_spin = handle_mo_item(start_mo)
                final_mo_num, final_spin = handle_mo_item(start_mo)
                last_exc_state = excited_states[-1]
                last_exc_state.add_mo_transition(start_mo_num, to_or_from,
                                                 final_mo_num, ci_coeff,
                                                 start_spin=start_spin,
                                                 final_spin=final_spin)
            # stop this matching if a blank line is encountered
            if line.strip() == "":
                matched_exc_state = False
            continue
        m_obj = re.match(EXC_LINE, line)
        if m_obj:
            excited_state = ExcitedState(*conv(m_obj.groups(), "issffff"))
            excited_state.mult = mult
            excited_states.append(excited_state)
            matched_exc_state = True

    return excited_states
