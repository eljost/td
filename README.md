# td.py

Parser for excited state calculations done with **Gaussian 09** (*td* keyword) or **TURBOMOLE** (only *escf* module).

By default the script expects **Gaussian 09** logs, but when the supplied filename is *escf.out* it will be parsed as **TURBOMOLE**-log.

The script uses a slightly modified version of Sergey Astanin's **tabulate** module (https://bitbucket.org/astanin/python-tabulate). Thanks to him.

## Installation
td.py requires:

    Python 3.x
    numpy
    
Optional packages are:

	python-docx
	argcomplete
	
The **argcomplete** module (https://pypi.python.org/pypi/argcomplete) has to be configured separately. When using **argcomplete** the script looks for files with *.out* and *.log* extensions.

## Usage
Display the help message with all available commands:

	td.py -h

Common usage examples are the generation of broadened spectra from calculated excitation energies and oscillator strenghts according to http://www.gaussian.com/g_whitepap/tn_uvvisplot.htm and this paper https://dx.doi.org/10.1002/chir.20733. Spectra can be generated in two ways: Unnormalized (ε in l mol⁻¹ cm⁻¹) or normalized with the brightest  peak set to 1 (in arbitrary units). The spectrum is printed to STDOUT.

### Spectrum generation

Normalized spectrum:

	./td.py [fn] --spectrum [from in nm] [to in nm]  > [outfn]
	
Unnormalized spectrum with --nnorm flag:

	./td.py [fn] --spectrum [from in nm] [to in nm] --nnorm > [outfn]
	
When used with the argument e2f the molecular extinction coefficients on the ordinate will be converted to a oscillator strength scale.

	./td.py [fn] --spectrum [from in nm] [to in nm] --e2f --nnorm > [outfn]

### Filtering

To investigate excited state optimizations it may be useful to split the output in chunks, where chunks should equal the number of calculated roots:

	./td.py [fn] --chunks [roots]

Only show transitions with an oscillator strength greater than or equal to a supplied threshold and sort by oscillator strength:
	
	./td.py [fn] --fthresh [thresh] --sf
