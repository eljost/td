# td.py

Parser for **td**-calculations done with **Gaussian 09**.

The script uses a slightly modified version of Sergey Astanin's **tabulate** module (https://bitbucket.org/astanin/python-tabulate). Thanks to him.

Installation
----------------
td.py requires:

    Python 3.x
    numpy
    
It also supports the **argcomplete** module (https://pypi.python.org/pypi/argcomplete),  but it's usage is optional. When using **argcomplete** the script looks for files with *.out* and *.log* extensions.

Usage
---------
Display the help message with all available commands:

	td.py -h

Common usage examples are the generation of broadened spectra from calculated excitation energies and oscillator strenghts according to http://www.gaussian.com/g_whitepap/tn_uvvisplot.htm and this paper https://dx.doi.org/10.1002/chir.20733. Spectra can be generated in two ways: Unnormalized (ε in l mol⁻¹ cm⁻¹) or normalized where the highest absorption equals 1 (in arbitrary units).

Normalized spectrum:

	./td.py [fn] --spectrum [from in nm] [to in nm]  > [outfn]
	
Unnormalized spectrum: Same as above, but with --nnorm flag.

	./td.py [fn] --spectrum [from in nm] [to in nm] --nnorm > [outfn]

To investigate excited state optimizations it may be useful to split the output in chunks, where chunks should equal the number of calculated roots:

	./td.py [fn] --chunks [roots]
