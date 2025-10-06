# ----- lyapunov.gnuplot -----
set terminal qt size 800,600 enhanced font "Arial,12"
set grid
set title "Evolucion del exponente de Lyapunov"
set xlabel "Tiempo (s)"
set ylabel "lambda(t)"
set key top right

plot "results/lyapunov.dat" using 1:4 every ::1 with lines lc rgb "blue" lw 2 title "lambda(t)"
pause -1 "Presiona Enter para salir"