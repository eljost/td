# td.py

Parser for excited state calculations done with **Gaussian 09** (*td* keyword) or **TURBOMOLE** (*escf* and *ricc2* module).

By default the script expects **Gaussian 09** logs, but when the supplied filename ends with *escf.out* or *ricc2* it will be parsed as **TURBOMOLE**-log of appropriate type.

The script uses a slightly modified version of Sergey Astanin's **tabulate** module (https://bitbucket.org/astanin/python-tabulate). Thanks to him.

## Installation
td.py requires:

    Python 3.x
    numpy
    simplejson
    
Optional packages are:

	python-docx
	argcomplete
	
The **argcomplete** module (https://pypi.python.org/pypi/argcomplete) has to be configured separately. When using **argcomplete** *td.py* looks for files with *.out* and *.log* extensions.

## Usage
Display the help message with all available commands:

	td.py -h

Common usage examples are the generation of broadened spectra from calculated excitation energies and oscillator strenghts according to http://www.gaussian.com/g_whitepap/tn_uvvisplot.htm and this paper https://dx.doi.org/10.1002/chir.20733. Spectra can be generated in two ways: Not normalized (ε in l mol⁻¹ cm⁻¹) or normalized with the brightest  peak set to 1 (ε/ε_max). The spectrum is printed to STDOUT.

### Verbose MO names 
Verbose MO names can optionally be loaded from a *mos.json* file with the following format:

	{
		"[MO number] [MO irrep] : "[verbose name]",
		e.g.,
		"61 a'" : "dπ₁",
		"39 a\"" : "dπ₂"
		
		! The last entry must NOT have a comma at the end
	}

### Spectrum generation

td.py exports the spectrum ε(x) for x in two different units: *nm* and *eV*. The first two blocks hold the spectrum and the oscillator strength impulses in *nm*, the third and fourth blocks hold the same data in *eV*. Within *gnuplot* the data blocks can be easily accessed by with the *index [id]* command.

#### Normalized spectrum with (ε/ε~max~) on the ordinate:

	./td.py [fn] --spectrum [from in nm] [to in nm]  > [outfn]
	
#### Non-normalized spectrum with ε on the ordinate:

	./td.py [fn] --spectrum [from in nm] [to in nm] --nnorm > [outfn]
	
#### Oscillator strength scale on the ordiante
When used with the argument  *\-\-e2f* the molecular extinction coefficients on the ordinate will be converted to an oscillator strength scale.

	./td.py [fn] --spectrum [from in nm] [to in nm] --e2f --nnorm > [outfn]

### Filtering

To investigate a **Gaussian 09** excited state optimization it may be useful to split the output in chunks, where chunks should equal the number of calculated roots at every step of the optimization:

	./td.py [fn] --chunks [roots]

Only show transitions with an oscillator strength greater than or equal to a supplied threshold and sort by oscillator strength:
	
	./td.py [fn] --fthresh [thresh] --sf

### Exporting
Several export-formats are available:

| Format | Argument | Comment |
| --------- | ------------- | ------------- |
| booktabs | \-\-booktabs | To be used in .tex-documents (doesn't export transitions and weights yet).|
| raw | \-\-raw | Export without any formatting. |
| docx | \-\-docx | Export to an .docx-document. Needs python-docx-module. |
| tiddlywiki | \-\-tiddly | Export as a tiddlywiki-table. |