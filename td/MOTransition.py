#!/usr/bin/env python3

class MOTransition:
    def __init__(self, start_mo, to_or_from, final_mo, ci_coeff,
                 contrib, start_spin, final_spin,
                 start_irrep, final_irrep):
        self.start_mo = int(start_mo)
        self.to_or_from = to_or_from
        self.final_mo = int(final_mo)
        self.ci_coeff = float(ci_coeff)

        self.contrib = float(contrib)
        self.start_spin = start_spin
        self.final_spin = final_spin
        self.start_irrep = start_irrep
        self.final_irrep = final_irrep

    def outstr(self):
        return "\t{0:>5d}{1} {2} {3} {4:>5d}{5} {6}" \
               "\t{7: 5.3f}\t{8:3.1%}".format(
                self.start_mo,
                self.start_irrep,
                self.start_spin[0],
                self.to_or_from,
                self.final_mo,
                self.final_irrep,
                self.final_spin[0],
                self.ci_coeff,
                self.contrib)

    def __str__(self):
        return self.outstr()

    def start_tpl(self):
        return (self.start_mo, self.start_irrep)

    def final_tpl(self):
        return (self.final_mo, self.final_irrep)
