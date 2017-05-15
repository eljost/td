#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

NM2EV = 1240.6691

class Spectrum:

    def __init__(self, excited_states):
        self.excited_states = excited_states

        wavelengths = [es.l for es in self.es]
        self.nm_range = np.array((min(wavelengths)-25,
                                  max(wavelengths)+100))
        """
        NM2EV = 1240.6691
        self.eV_range = NM2EV / self.nm_range
        """

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
        in_nm = np.stack((x, spectrum), axis=-1)
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

    def plot_eV(self):
        in_eV, osc_eV = self.eV
        print(osc_eV)
        fig, ax1 = plt.subplots()

        ax1.plot(in_eV[:,0], in_eV[:,1])
        from_x, to_x = ax1.get_xlim()
        ax1.set_xlim(to_x, from_x)
        ax1.set_xlabel("E / eV")
        ax1.set_ylabel("ε / mol cm⁻¹ l⁻¹")

        ax2 = ax1.twinx()
        ax2.stem(osc_eV[:,0], osc_eV[:,1], markerfmt=" ", basefmt=" ")

        plt.show()
