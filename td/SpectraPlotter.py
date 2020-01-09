#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from td.peakdetect import peakdetect
from td.constants import NM2EV


class SpectraPlotter:

    def __init__(self, spectra, unit, peaks=False, enum=False):
        self.spectra = spectra
        self.unit = unit
        self.peaks = peaks
        self.enum = enum

        self.broadened = [getattr(spectrum, self.unit) for spectrum
                          in self.spectra]
        self.fig, self.ax = plt.subplots()
        self.ax2 = self.ax.twinx()
        self.colors = mpl.rcParams["axes.prop_cycle"]

    def plot_spectrum(self, spectrum, color):
        broad_spec, osc_spec = getattr(spectrum, self.unit)

        self.ax.plot(broad_spec[:,0], broad_spec[:,1],
                     label=f"{spectrum.name}", **color)

        # Enumerate stick plot
        if self.enum:
            for i, stick in enumerate(osc_spec, 1):
                x, y = stick 
                self.ax2.text(x, y*1.1, f"$S_{i}$")

        if self.peaks:
            peak_inds = spectrum.get_peak_inds(broad_spec)
            peaks = broad_spec[peak_inds]
            print(f"Peaks from {spectrum.name}:")
            for i, peak in enumerate(peaks):
                energy, epsilon = peak[:2]
                print("{:2d}:\tλ={:5.2f} {}, ε={:8.0f}".format(i, energy,
                                                             self.unit, epsilon))
                xytext = peak[:2] * (1, 1.05)
                self.ax.annotate("{}".format(i), xy=peak[:2], xytext=xytext,
                                 horizontalalignment="center")
            self.ax.plot(peaks[:,0], peaks[:,1], "o", **color)

        stemlines =  self.ax2.stem(osc_spec[:,0], osc_spec[:,1],
                                   markerfmt=" ", basefmt=" ")
        plt.setp(stemlines, 'color', color["color"])


    def plot(self):

        xlabel = "E / {}".format(self.unit)

        for spectrum, color in zip(self.spectra, self.colors):
            self.plot_spectrum(spectrum, color)

        if self.unit == "eV":
            from_x, to_x = self.ax.get_xlim()
            self.ax.set_xlim(to_x, from_x)
        elif self.unit == "nm":
            self.ax.set_xlim(self.spectra[0].nm_range)

        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel("ε / mol cm⁻¹ l⁻¹")
        self.ax2.set_ylabel("f")

        if len(self.spectra) == 1:
            self.fig.suptitle(self.spectra[0].name)
        else:
            self.ax.legend()

        #from_y2, to_y2 = ax2.get_ylim()
        #to_y2 = max(to_y2, 0.5)
        #ax2.set_ylim(from_y2, to_y2)

        plt.show()
