#!/usr/bin/env python3

import logging
import math
import re

from td.helper_funcs import IRREPS_REPL
from td.MOTransition import MOTransition

class ExcitedState:
    def __init__(self, id, spin, spat, dE, l, f, s2):
        self.id = int(id)
        self.id_sorted = self.id
        self.spin = spin
        self.spat = spat
        self.dE = dE
        self.l = l
        self.f = f
        self.s2 = s2
        self.rr_weight = None

        self.mo_transitions = list()
        self.irrep = self.spat
        self.normalize_irrep()

    def normalize_irrep(self):
        for key in IRREPS_REPL:
            self.irrep = re.sub(key, IRREPS_REPL[key], self.irrep)

    def as_list(self, attrs=None):
        if not attrs:
            attrs = ("id",
                     "id_sorted",
                     "spin",
                     "spat",
                     "dE",
                     "l",
                     "f",
                     "s2")
        return [getattr(self, a) for a in attrs]

    def add_mo_transition(self, start_mo, to_or_from, final_mo, ci_coeff,
                          start_spin="α", final_spin="α",
                          contrib=0.0,
                          start_irrep="a", final_irrep="a"):
        start_spin = "α" if (start_spin == "a") else "β"
        final_spin = "α" if (final_spin == "a") else "β"

        self.mo_transitions.append(MOTransition(
            start_mo, to_or_from, final_mo, ci_coeff,
            contrib, start_spin, final_spin, start_irrep, final_irrep)
        )

    def get_start_mos(self):
        return set([mo_tr.start_mo for mo_tr in self.mo_transitions])

    def get_final_mos(self):
        return set([mo_tr.final_mo for mo_tr in self.mo_transitions])

    def has_mo_transition(self, start_mo, final_mo):
        return bool(
            [mo_exc for mo_exc in self.mo_transitions
             if (start_mo == mo_exc.start_mo) and
                (final_mo == mo_exc.final_mo)]
        )

    def is_singlet(self):
        return self.mult == 1

    def calculate_contributions(self):
        for mo_trans in self.mo_transitions:
            if (mo_trans.contrib is not None and
                not math.isclose(mo_trans.contrib, 0.0)):
                continue
            contrib = mo_trans.ci_coeff**2
            if self.is_singlet():
                contrib *= 2
            mo_trans.contrib = contrib

    def correct_backexcitations(self):
        # Check if there are any back-excitations, e.g.
        # 89B <- 90B
        back_transitions = [bt for bt in self.mo_transitions
                            if bt.to_or_from == "<-"]
        for bt in back_transitions:
            # Find corresponding transition from final_mo -> start_mo
            trans_to_correct = [t for t in self.mo_transitions
                                if (t.start_mo == bt.start_mo and
                                    t.final_mo == bt.final_mo and
                                    t.to_or_from == "->")]
            if trans_to_correct:
                assert(len(trans_to_correct) == 1)
                trans_to_correct = trans_to_correct[0]
                # Correct contribution of trans_to_correct

                trans_to_correct.contrib -= bt.contrib

    def print_mo_transitions(self, verbose_mos):
        for mot in self.mo_transitions:
            # Suppresss backexcitations like
            # 89B <- 90B
            if mot.to_or_from == "<-":
                continue
            print(mot.outstr())
            if verbose_mos:
                try:
                    start_mo_verbose = verbose_mos[mot.start_tpl()]
                    final_mo_verbose = verbose_mos[mot.final_tpl()]
                    print("\t\t{0} -> {1}".format(
                        start_mo_verbose,
                        final_mo_verbose
                    ))
                except KeyError as err:
                    logging.warning("Verbose MO name for {} {}"
                                    " missing!".format(*err.args[0]))

    def suppress_low_ci_coeffs(self, thresh):
        # Check if excitation lies below ci coefficient threshold
        # if it does don't add this excitation
        contrib_thresh = thresh**2 * 2
        to_del = list()
        for mo_trans in self.mo_transitions:
            if mo_trans.contrib:
                if mo_trans.contrib <= contrib_thresh:
                    to_del.append(mo_trans)
                continue
            if abs(mo_trans.ci_coeff) <= thresh:
                to_del.append(mo_trans)
        for td in to_del:
            self.mo_transitions.remove(td)

    def calc_rr_weight(self, rr_exc):
        l_in_cm = 10**7 / self.l
        rr_ex_in_cm = 10**7 / rr_exc
        G = 1500j
        self.rr_weight = self.f * abs(G/(l_in_cm-rr_ex_in_cm-G))

    def update_irreps(self):
        """Find out which irreps are important for the MOTransitions."""
        start_irreps = [mo_trans.start_irrep for mo_trans
                        in self.mo_transitions]
        final_irreps = [mo_trans.final_irrep for mo_trans
                        in self.mo_transitions]
        self.irreps = set(start_irreps + final_irreps)

        self.mo_trans_per_irrep = dict()
        for irrep in self.irreps:
            start_mos = [mot.start_mo for mot in self.mo_transitions
                         if mot.start_irrep == irrep]
            final_mos = [mot.final_mo for mot in self.mo_transitions
                         if mot.final_irrep == irrep]
            unique_mos = set(start_mos + final_mos)
            self.mo_trans_per_irrep[irrep] = unique_mos

    """
    # find lowest orbital from where an excitation originates
    min_mo = min(itertools.chain(*[exc_state.get_start_mos()
                                   for exc_state in excited_states]))
    # find highest lying orbital where an excitation ends
    max_mo = max(itertools.chain(*[exc_state.get_final_mos()
                                   for exc_state in excited_states]))
    """

    def __str__(self):
        return "#{0} {1} eV f={2}".format(self.id, self.dE, self.f)
