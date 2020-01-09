#!/usr/bin/env gnuplot

set terminal pdfcairo enhanced linewidth 2
set output "nm.pdf"

unset key

abs_ylabel = "{/Symbol e} l mol⁻¹ cm⁻¹"
abs_norm_ylabel = "{/Symbol e} / {/Symbol e}_{max}"

set ytics nomirror
set y2tics
set y2range [0:0.75]
set y2label "f"
set ylabel abs_ylabel
set xlabel "E / nm"
set yrange [0:]

set style line 1 lc rgb "#0571b0"

plot "nm.spec" w l ls 1, \
 "osc_nm.spec" with impulses axes x1y2 notitle ls 1

set output "nm_norm.pdf"
set ylabel abs_norm_ylabel
plot "nm.spec" u 1:3 w l ls 1, \
 "osc_nm.spec" with impulses axes x1y2 notitle ls 1 #, \
 "< awk '(NR>2){print;}' path" u 1:5 w l axes x1y1

set output "eV.pdf"

set xrange [] reverse

set xlabel "E / eV"
set ylabel abs_ylabel
plot "eV.spec" w l ls 1, \
 "osc_eV.spec" with impulses axes x1y2 notitle ls 1

set output "eV_norm.pdf"
set ylabel abs_norm_ylabel
plot "eV.spec" u 1:3 w l ls 1, \
 "osc_eV.spec" with impulses axes x1y2 notitle ls 1
