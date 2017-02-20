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
    lines = text.split("\n")
    excited_states = list()
    involved_mos = dict()

    matched_exc_state = False
    for line in lines:
        line = line.strip()
        # look for initial and final MO and corresponding CI-coefficient
        if matched_exc_state:
            match_obj = re.match(TRS_LINE, line)
            if match_obj:
                group_list = list(match_obj.groups())
                conv_mo_excitation = conv(group_list, "sssf")
                last_exc_state = excited_states[-1]
                last_exc_state.add_mo_transition(*conv_mo_excitation)
                try:
                    involved_mos[last_exc_state.id].append(
                            conv_mo_excitation)
                except KeyError:
                    involved_mos[last_exc_state.id] = [conv_mo_excitation, ]
            # stop this matching if a blank line is encountered
            if line.strip() == "":
                matched_exc_state = False
            continue
        m_obj = re.match(EXC_LINE, line)
        if m_obj:
            excited_state = ExcitedState(*conv(m_obj.groups(), "issffff"))
            excited_states.append(excited_state)
            matched_exc_state = True

    return excited_states, involved_mos


if __name__ == "__main__":
    fn = "/home/carpx/Dokumente/projekte/runobpb/runobpb/01_runobpb_cam_cpcm_td.log"
    with open(fn) as handle:
        text = handle.read()
    es, mos = parse_tddft(text)
    print(es)
