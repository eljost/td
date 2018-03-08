#!/usr/bin/env python3

"""This script plots two gradually mixed spectra. It can be used
to visualize the spectra changes over the course of a chemical reaction,
e.g. a photoreaction."""

import argparse
import sys

import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("spec1",
                        help="First spectrum, produced by td.py --savenm")
    parser.add_argument("spec2",
                        help="Second spectrum, produced by td.py --savenm")
    parser.add_argument("--xlim", nargs=2, type=int,
                        help="Range in nm to use for plotting.")
    parser.add_argument("--ylim", nargs=2, type=int,
                        help="Range for the ordinate.")
    parser.add_argument("--legend", action="store_true")

    return parser.parse_args(args)


def mix_vectors(vec1, vec2, fact):
    return fact*vec1+(1-fact)*vec2


def assert_nm_range(nm, step=0.5):
    diff = np.ediff1d(nm)
    np.testing.assert_allclose(diff, np.full_like(diff, step),
        err_msg="Spectrum has to be given in 0.5 nm steps!")


def run():
    args = parse_args(sys.argv[1:])
    spec1 = np.loadtxt(args.spec1)
    spec2 = np.loadtxt(args.spec2)
    nm1 = spec1[:,0]
    nm2 = spec2[:,0]

    assert_nm_range(nm1)
    assert_nm_range(nm2)

    nm_min = min(nm1.min(), nm2.min())
    nm_max = max(nm1.max(), nm2.max())
    # Add 0.5 because the last point is exclusive
    nm_new = np.arange(nm_min, nm_max+0.5, 0.5)
    # Find indices of nm1 and nm2 in nm_new
    nm1_start = np.where(nm_new == nm1[0])[0][0]
    nm2_start = np.where(nm_new == nm2[0])[0][0]
    # Exctinction coeffs. of the convoluted spectra
    eps1 = spec1[:,1]
    eps2 = spec2[:,1]
    spec1_full = np.zeros_like(nm_new)
    spec2_full = np.zeros_like(nm_new)
    spec1_full[nm1_start:nm1_start+eps1.size] = eps1
    spec2_full[nm2_start:nm2_start+eps2.size] = eps2

    fig, ax = plt.subplots()

    frames = 21
    colors = plt.cm.coolwarm(np.linspace(0.1,0.9,frames))

    factors = np.linspace(1, 0, frames, endpoint=True)
    for i in range(frames):
        fact = factors[i]
        color = colors[i]
        mixed = mix_vectors(spec1_full, spec2_full, fact)
        label = None
        if i % 2 == 0:
            label = f"{1-fact:.0%} Product"
        ax.plot(nm_new, mixed, c=color, label=label)
    if args.xlim:
        ax.set_xlim(*args.xlim)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    if args.legend:
        ax.legend()
    ax.set_xlabel("λ / nm")
    ax.set_ylabel("ε / l mol⁻¹ cm⁻¹")
    plt.show()


if __name__ == "__main__":
    run()
