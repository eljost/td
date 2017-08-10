#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

from td.helper_funcs import EV2NM, FLOAT_RE
from td.ExcitedState import ExcitedState


def parse_tddft(text):
    # Use transition electric dipole moments
    abs_re = "VIA TRANSITION ELECTRIC DIPOLE MOMENTS(.+?)-\s*A"
    abs_text = re.search(abs_re, text, re.DOTALL).groups()[0]
    line_re = ("(\d+)", FLOAT_RE, FLOAT_RE, FLOAT_RE, ".+\n")
    line_re = "\s*".join(line_re)
    states = re.findall(line_re, abs_text)
    states = [(int(id), float(l), float(f)) for id, ecm, l, f in states]

    states_moc_re = "STATE\s*\d+:\s*E=(.+?)\n\n"
    states_moc = re.findall(states_moc_re, text, re.DOTALL)
    assert(len(states) == len(states_moc))

    excited_states = list()
    moc_re = ("(\d+)(a|b)", "->", "(\d+)(a|b)", ":", FLOAT_RE, "\(c=",
              FLOAT_RE)
    moc_re = "\s*".join(moc_re)
    # No spatial symmetry in ORCA
    sym = "a"
    start_irrep = "a"
    final_irrep = "a"
    spin = "???"
    to_or_from = "->"
    for (id, l, f), state_moc in zip(states, states_moc):
        ee = EV2NM / l
        exc_state = ExcitedState(id, spin, sym, ee, l, f, "???")
        excited_states.append(exc_state)
        mocs = re.findall(moc_re, state_moc)
        for (start_mo, start_spin, final_mo, final_spin,
             percent, coeff) in mocs:
            start_mo = int(start_mo)
            final_mo = int(final_mo)
            percent = float(percent)
            coeff = float(coeff)
            exc_state.add_mo_transition(start_mo,
                                        to_or_from,
                                        final_mo,
                                        ci_coeff=coeff,
                                        contrib=percent,
                                        start_spin=start_spin,
                                        final_spin=final_spin,
                                        start_irrep=start_irrep,
                                        final_irrep=final_irrep)

    return excited_states


def parse_nto_block(text):
    state_re = "STATE\s*(\d+)"
    state = int(re.search(state_re, text).groups()[0])
    nto_contrib_re = "(\d+)([ab])\s*->\s*(\d+)([ab])\s*:\s*n=\s*([\d\.]+)"
    nto_contribs = re.findall(nto_contrib_re, text)
    nto_contribs = [(int(from_nto), from_spin,
                     int(to_nto), to_spin,
                     float(nto_weight))
                    for from_nto, from_spin, to_nto, to_spin, nto_weight
                    in nto_contribs]
    return state, nto_contribs


def parse_ntos(text):
    nto_block_re = "(NATURAL TRANSITION ORBITALS FOR STATE.+?-+.+?-----)"
    nto_blocks = re.findall(nto_block_re, text, re.DOTALL)
    parsed_nto_blocks = [parse_nto_block(ntob) for ntob in nto_blocks]
    return parsed_nto_blocks

if __name__ == "__main__":
    fn = "11_b3lyp35_cpcm_120_tddft.out"
    fn = "/scratch/nbdcorm/koligand/11_b3lyp35_cpcm_120_tddft/11_b3lyp35_cpcm_120_tddft.out"
    with open(fn) as handle:
        text = handle.read()
    parsed_nto_blocks = parse_ntos(text)
    states, nto_contribs = zip(*parsed_nto_blocks)
    print(states[0])
    print(nto_contribs[0])
