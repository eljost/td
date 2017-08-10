#!/usr/bin/env gnuplot

set terminal pdfcairo enhanced linewidth 2
set output "eV.pdf"

set ytics nomirror
set y2tics
set y2range [0:1]
set y2label "f"
set ylabel "{/Symbol e} l mol⁻¹ cm⁻¹"
set xlabel "E / eV"

set style line 1 lc rgb "#0571b0"

plot "eV.spec" w l ls 1, \
 "osc_eV.spec" with impulses axes x1y2 notitle ls 1
