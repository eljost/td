#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Johannes Steinmetzer, 2017
# PYTHON_ARGCOMPLETE_OK

import argparse
import itertools
import logging
import os
import re
import shutil
import sys

import numpy as np
import simplejson as json

from td.helper_funcs import chunks, THIS_DIR
from td.ExcitedState import ExcitedState
from td.export import *
import td.parser.gaussian as gaussian
import td.parser.orca as orca
import td.parser.turbomole as turbo
from td.Spectrum import Spectrum

# Optional modules
try:
    import argcomplete
except ImportError:
    pass


def is_orca(text):
    orca_re = "\* O   R   C   A \*"
    return re.search(orca_re, text)


def set_ntos(excited_states, ntos):
    for es, (state, nto_contribs) in zip(excited_states, ntos):
        assert(es.id == state)
        es.mo_transitions = list()
        for from_nto, from_spin, to_nto, to_spin, nto_weight in nto_contribs:
            es.add_mo_transition(start_mo=from_nto,
                                 to_or_from="->",
                                 final_mo=to_nto,
                                 ci_coeff=0,
                                 start_spin=from_spin,
                                 final_spin=to_spin,
                                 contrib=nto_weight)
    return excited_states


def determine_program(args):
    fn = args.file_name
    with open(fn) as handle:
        text = handle.read()

    # Parse outputs
    if fn.endswith("escf.out"):
        # TURBOMOLE escf
        excited_states = turbo.parse_escf(text)
    elif fn.endswith("ricc2.out"):
        # TURBOMOLE ricc2
        excited_states = turbo.parse_ricc2(text)
    elif is_orca(text):
        excited_states = orca.parse_tddft(text)
        if args.ntos:
            ntos = orca.parse_ntos(text)
            excited_states = set_ntos(excited_states, ntos)
    else:
        # Gaussian
        excited_states = gaussian.parse_tddft(text)

    return excited_states


def logs_completer(prefix, **kwargs):
    return [path for path in os.listdir(".") if path.endswith(".out") or
            path.endswith(".log")]


def parse_args(args):
    parser = argparse.ArgumentParser(
            "Displays output from Gaussian-td-calculation,"
            " sorted by oscillator strength f."
    )

    parser.add_argument("--show", metavar="n", type=int,
                        help="Show only the first n matching excitations.")
    parser.add_argument("--only-first", metavar="only_first", type=int,
                        help="Only consider first n excited states.")
    parser.add_argument("--range", metavar="start_end", nargs="+", type=float,
                        help="Show only excited states accessible in this "
                             "wavelength range (e.g. 400 450).")

    sorting_group = parser.add_mutually_exclusive_group()
    sorting_group.add_argument("--sf", action="store_true",
                               help="Sort by oscillator strength.")
    sorting_group.add_argument("--se", action="store_true",
                               help="Sort by energy.")

    parser.add_argument("--start-mos", dest="start_mos", type=str, nargs="+",
                        help="Show only transitions from this MO.")
    parser.add_argument("--final-mos", dest="final_mos", type=str, nargs="+",
                        help="Show only transitions to this MO.")
    parser.add_argument("--start-final-mos", dest="start_final_mos",
                        type=str, nargs="+", help="(Number of) MO pair(s). "
                        "Only transitions from [start mo] to [final mo] "
                        "are shown.")
    parser.add_argument("--raw", action="store_true",
                        help="Just print the data, without the table.")
    parser.add_argument("--by-id", dest="by_id", type=int,
                        help="Display excited state with specific id.")
    parser.add_argument("--summary", action="store_true",
                        help="Print summary to stdout.")
    parser.add_argument("--ci-coeff", dest="ci_coeff", type=float,
                        default=0.2, help="Only consider ci coefficients "
                        "not less than.")
    parser.add_argument("--spectrum", dest="spectrum", action="store_true",
                        help="Calculate the UV spectrum from the TD "
                        "calculation (FWHM = 0.4 eV).")
    parser.add_argument("--e2f", dest="e2f", action="store_true",
                        help="Used with spectrum. Converts the molecular "
                        "extinctions coefficients on the ordinate to "
                        "oscillator strengths.")
    """
    parser.add_argument("--hi", dest="highlight_impulses", type=int,
                        nargs="+", help="List of excitations. Their "
                        "oscillator strength bars will be printed separatly.")
    """
    parser.add_argument("--norm", type=int,
                        help="Normalize the calculated spectrum.")
    parser.add_argument("--irrep", dest="irrep",
                        help="Filter for specific irrep.")
    parser.add_argument("--booktabs", dest="booktabs", action="store_true",
                        help="Output table formatted for use with the latex-"
                        "package booktabs.")
    parser.add_argument("--rrexc", type=float,
                        help="Excitation wavelength for resonance raman.")
    parser.add_argument("--rrthresh", type=float, default=1e-2,
                        help="Threshold for RR weight.")
    parser.add_argument("--fthresh", type=float, default=0.0,
                        help="Only show transitions with oscillator strengths"
                        " greater than or equal to the supplied threshold.""")
    parser.add_argument("--chunks", type=int, default=0,
                        help="Split the output in chunks. Useful for "
                        "investigating excited state optimizations. Don't use "
                        "with --sf or --raw.")
    parser.add_argument("--docx", action="store_true",
                        help="Output the parsed data as a table into a "
                        ".docx document.")
    parser.add_argument("--tiddly", action="store_true",
                        help="Output the parsed data in Tiddlywiki-table"
                        "format.")
    parser.add_argument("--theodore", action="store_true",
                        help="Output HTML with NTO-picture from THEOdore.")
    parser.add_argument("--nosym", action="store_true",
                        help="Assign all excited states to the 'a' irrep.")
    parser.add_argument("--plot", choices=["eV", "nm"],
                        help="Plot the spectrum with matplotlib.")
    parser.add_argument("--peaks", action="store_true", default=False,
                        help="Detect peaks.")
    parser.add_argument("--ntos", action="store_true", default=False,
                        help="Parse NTOs from an ORCA log.")

    # Use the argcomplete module for autocompletion if it's available
    if "argcomplete" in sys.modules:
        parser.add_argument("file_name", metavar="fn",
                            help="File to parse.").completer = logs_completer
        argcomplete.autocomplete(parser)
    else:
        parser.add_argument("file_name", metavar="fn",
                            help="File to parse.")
    return parser.parse_args(args)

def run():
    args = parse_args(sys.argv[1:])

    logging.info("Only considering transitions  with "
                 "CI-coefficients >= {}:".format(args.ci_coeff))

    try:
        mo_data_fn = "mos.json"
        with open(mo_data_fn) as handle:
            json_data = json.load(handle)

        verbose_mos = dict()
        for key in json_data:
            mo, irrep = key.split()
            #mo = int(mo)
            verbose_mos[(mo, irrep)] = json_data[key]
    except IOError:
        logging.warning("Couldn't find verbose MO-names"
                        " in mos.json")
        verbose_mos = None


    excited_states = determine_program(args)
    
    for exc_state in excited_states:
        exc_state.calculate_contributions()
        exc_state.correct_backexcitations()
        exc_state.suppress_low_ci_coeffs(args.ci_coeff)
        exc_state.update_irreps()

    spectrum = Spectrum(excited_states)

    if args.nosym:
        for es in excited_states:
            es.spat = "a"
            es.irrep = "a"

    if args.only_first:
        excited_states = excited_states[:args.only_first]

    if args.plot == "eV":
        spectrum.plot_eV(title=args.file_name, with_peaks=args.peaks)
    elif args.plot == "nm":
        spectrum.plot_nm(title=args.file_name, with_peaks=args.peaks)

    if args.by_id:
        try:
            exc_state = excited_states[args.by_id-1]
            print_table([exc_state, ])
            exc_state.print_mo_transitions(verbose_mos)
        except IndexError:
            print("Excited state with id #{} not found.".format(args.by_id))
        sys.exit()

    """
    !
    ! DO FILTERING/SORTING HERE
    !
    """

    if args.irrep:
        excited_states = [es for es in excited_states if es.spat == args.irrep]

    if args.start_mos:
        states = set()
        for start_mo in args.start_mos:
            states.update([exc_state for exc_state in excited_states
                           if start_mo in exc_state.get_start_mos()])
        excited_states = states
    if args.final_mos:
        states = set()
        for final_mo in args.final_mos:
            states.update([exc_state for exc_state in excited_states
                           if final_mo in exc_state.get_final_mos()])
        excited_states = states
    if args.start_final_mos:
        sf_mos = args.start_final_mos
        if (len(sf_mos) % 2) != 0:
            sys.exit("Need an even number of arguments for "
                     "--start-final-mos, not an odd number.")
        states = set()
        pairs = [(sf_mos[i], sf_mos[i+1]) for i in range(len(sf_mos) / 2)]
        for start_mo, final_mo in pairs:
            states.update(
                [exc_state for exc_state in excited_states if
                 exc_state.has_mo_transition(start_mo, final_mo)]
            )
        excited_states = states

    # Sort by oscillator strength if requested.
    if args.sf:
        excited_states = sorted(excited_states,
                                key=lambda exc_state: -exc_state.f)
    # Sort by energy if requested.
    if args.se:
        excited_states = sorted(excited_states,
                                key=lambda es: -es.l)
    for i, es in enumerate(excited_states, 1):
        es.id_sorted = i

    # Only show excitations in specified wavelength-range
    if args.range:
        # Only lower threshold specified (energy wise)
        if len(args.range) is 1:
            end = args.range[0]
            excited_states = [exc_state for exc_state in excited_states
                              if (exc_state.l >= end)]
        elif len(args.range) is 2:
            start, end = args.range
            excited_states = [exc_state for exc_state in excited_states
                              if (start <= exc_state.l <= end)]
        else:
            raise Exception("Only 1 or 2 arguments allowed for --range!")

    if args.fthresh:
        excited_states = [exc_state for exc_state in excited_states
                          if (exc_state.f >= args.fthresh)]

    excited_states = list(excited_states)[:args.show]

    """
    !
    ! DON'T DO FILTERING/SORTING BELOW THIS LINE
    !
    """

    # Find important irreps
    irreps = set(itertools.chain(*[exc_state.irreps for exc_state in
                                   excited_states]))
    min_max_mos = dict()
    for irrep in irreps:
        min_max_mos[irrep] = list()
    for exc_state in excited_states:
        for irrep in exc_state.irreps:
            min_max_mos[irrep].extend(exc_state.mo_trans_per_irrep[irrep])
    for irrep in min_max_mos:
        mos = set(min_max_mos[irrep])
        min_max_mos[irrep] = (min(mos), max(mos))

    # Convert remaining excitations to a list so it can be printed by
    # the tabulate module
    as_list = [exc_state.as_list() for exc_state in excited_states]

    """
    !
    ! PRINTING BELOW THIS LINE
    !
    """

    if args.docx:
        as_docx(excited_states, verbose_mos)
    if args.tiddly:
        as_tiddly_table(excited_states, verbose_mos)
    if args.theodore:
        as_theodore(excited_states, args.file_name)
    if args.spectrum:
        # Starting and ending wavelength of the spectrum to be calculated
        in_nm, osc_nm = spectrum.nm
        in_eV, osc_eV = spectrum.eV
        if args.norm or (args.norm == 0):
            peak_inds = spectrum.get_peak_inds(in_nm)[args.norm]
            nm_peaks = in_nm[peak_inds]
            eV_peaks = in_eV[peak_inds]
            in_nm[:,2] = in_nm[:,1] / nm_peaks[1]
            in_eV[:,2] = in_eV[:,1] / eV_peaks[1]
        out_fns = ["nm.spec", "osc_nm.spec", "eV.spec", "osc_eV.spec"]
        for out_fn, spec in zip(out_fns, (in_nm, osc_nm, in_eV, osc_eV)):
            np.savetxt(out_fn, spec)
        gnuplot_tpl = os.path.join(THIS_DIR, "templates", "gnuplot.plt")
        shutil.copy(gnuplot_tpl, "gnuplot.plt")

    # Dont print the pretty table when raw output is requested
    # Don't print anything after the summary
    if args.raw:
        for exc_state in as_list:
            print("\t".join([str(item) for item in exc_state]))
    elif args.chunks > 0:
        for i, chunk in enumerate(chunks(excited_states, args.chunks), 1):
            print("### Chunk {} ###".format(i))
            print_table(chunk)
            print()
    elif args.summary:
        for exc_state in excited_states:
            print_table([exc_state, ])
            exc_state.print_mo_transitions(verbose_mos)
            print("")
    elif args.booktabs:
        as_booktabs(excited_states)
    else:
        # Print as pretty table with header information
        print_table(excited_states)

    if args.rrexc:
        for es in excited_states:
            es.calc_rr_weight(args.exc)
        rr_weights = [(es.id, es.rr_weight) for es in excited_states
                      if es.rr_weight >= args.rrthresh]
        print(tabulate(rr_weights))

    for irrep in irreps:
        min_mo, max_mo = min_max_mos[irrep]
        print("Irrep {}: MOs {} - {}".format(irrep,
                                             min_mo,
                                             max_mo))

if __name__ == "__main__":
    run()
