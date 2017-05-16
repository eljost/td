#!/usr/bin/env gnuplot

set terminal pdfcairo enhanced linewidth 2
set output "nm.pdf"

unset key

set ytics nomirror
set y2tics
set y2range [0:1]
set y2label "f"
set ylabel "{/Symbol e} l mol⁻¹ cm⁻¹"
set xlabel "E / nm"

set style line 1 lc rgb "#0571b0"

plot "nm.spec" w l ls 1, \
 "osc_nm.spec" with impulses axes x1y2 notitle ls 1

set output "eV.pdf"

set xlabel "E / eV"
plot "eV.spec" w l ls 1, \
 "osc_eV.spec" with impulses axes x1y2 notitle ls 1
