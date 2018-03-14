#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from td.peakdetect import peakdetect

NM2EV = 1240.6691

class Spectrum:

    def __init__(self, name, excited_states, gs_energy=None):
        self.name = name
        self.excited_states = excited_states

        wavelengths = [es.l for es in self.es]
        self.nm_range = np.array((int(min(wavelengths))-25,
                                  int(max(wavelengths))+100))
        """
        NM2EV = 1240.6691
        self.eV_range = NM2EV / self.nm_range
        """
        self.gs_energy = gs_energy

    @property
    def es(self):
        return self.excited_states

    def gauss_uv_band(self, x, osc, x_i):
        return (1.3062974e8 * osc / (1e7 / 3099.6) *
                np.exp(-((1. / x - 1. / x_i) / (1. / 3099.6))**2))

    @property
    def nm(self):
        return self.convolute(*self.nm_range)

    @property
    def eV(self):
        in_nm, osc_nm = self.nm
        in_eV = in_nm
        in_eV[:,0] = NM2EV / in_nm[:,0]
        osc_in_eV = osc_nm 
        osc_in_eV[:,0] = NM2EV / osc_in_eV[:,0]
        return in_eV, osc_in_eV

    def convolute(self, from_nm, to_nm):
        # According to:
        # http://www.gaussian.com/g_whitepap/tn_uvvisplot.htm
        # wave lengths and oscillator strengths
        # E(eV) = 1240.6691 eV * nm / l(nm)
        NM2EV = 1240.6691

        osc_nm = np.array([(es.l, es.f) for es in self.excited_states])
        x = np.arange(from_nm, to_nm, 0.5)
        spectrum = list()
        for l in x:
            spectrum.append(np.sum([self.gauss_uv_band(l, osc, l_i)
                            for l_i, osc in osc_nm]))
        spectrum = np.array(spectrum)
        spectrum_norm = spectrum / spectrum.max()
        in_nm = np.stack((x, spectrum, spectrum_norm), axis=-1)
        return in_nm, osc_nm

        """
        if e2f:
            spectrum /= 40490.05867167
        if not nnorm:
            spectrum = spectrum / spectrum.max()
        """

        """
        Used for printing also the gauss bands of the
        n-highest transitions
        # Sort by f
        fli_sorted = sorted(fli, key=lambda tpl: -tpl[0])
        highest_fs = fli_sorted[:15]
        print
        x = np.arange(200, 600, 0.5)
        highest_bands = list()
        for f, l_i in highest_fs:
            band = lambda l: gauss_uv_band(l, f, l_i)
            calced_band = band(x)
            calced_band = band(x) / calced_band.max() * f * 3
            highest_bands.append(calced_band)
            #for xi, ai in zip(x, calced_band):
            #    print xi, ai
            #print
        highest_bands_headers = tuple(['"l_i = {} nm, f = {}"'.format(
            l_i, f) for f, l_i in highest_fs])
        headers = ('"l in nm"', "Sum") + highest_bands_headers
        wargel = zip(x, spectrum, *highest_bands)
        print tabulate(wargel, headers=headers, tablefmt="plain")
        """

    def plot_eV(self, title="", with_peaks=False):
        in_eV, osc_eV = self.eV
        unit = "eV"
        self.plot(in_eV, osc_eV, unit, title=title, reverse_x=True,
                  with_peaks=with_peaks)

    def plot_nm(self, title="", with_peaks=False):
        in_nm, osc_nm = self.nm
        unit = "nm"
        self.plot(in_nm, osc_nm, unit, title=title, with_peaks=with_peaks)

    def plot(self, conv_spectrum, osc, unit, title="", reverse_x=False,
             with_peaks=None):
        fig, ax1 = plt.subplots()
        fig.suptitle(title)
        xlabel = "E / {}".format(unit)

        ax1.plot(conv_spectrum[:,0], conv_spectrum[:,1])

        if with_peaks:
            peak_inds = self.get_peak_inds(conv_spectrum)
            peaks = conv_spectrum[peak_inds]
            for i, peak in enumerate(peaks):
                energy, epsilon = peak[:2]
                print("{:2d}:\tλ={:5.0f} {}, ε={:8.0f}".format(i, energy,
                                                             unit, epsilon))
                xytext = peak[:2] * (1, 1.05)
                ax1.annotate("{}".format(i), xy=peak[:2], xytext=xytext,
                             horizontalalignment="center")
            ax1.plot(peaks[:,0], peaks[:,1], "ro")

        if reverse_x:
            from_x, to_x = ax1.get_xlim()
            ax1.set_xlim(to_x, from_x)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel("ε / mol cm⁻¹ l⁻¹")

        ax2 = ax1.twinx()
        ax2.stem(osc[:,0], osc[:,1], markerfmt=" ", basefmt=" ")

        from_y2, to_y2 = ax2.get_ylim()
        to_y2 = max(to_y2, 0.5)
        ax2.set_ylim(from_y2, to_y2)
        ax2.set_ylabel("f")

        plt.show()

    def get_peak_inds(self, conv_spectrum, lookahead=25):
        conv_spectrum_ys = conv_spectrum[:,1]
        max_peaks, min_peaks = peakdetect(conv_spectrum_ys, lookahead=lookahead)
        return np.array(max_peaks)[:,0].astype(int)


    def write_nm(self):
        in_nm, _ = self.nm
        nm_name = f"{self.name}_nm.dat"
        np.savetxt(nm_name, in_nm)


    def __str__(self):
        return f"Spectrum({len(self.excited_states)} states)"

