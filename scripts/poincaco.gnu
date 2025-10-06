set encoding utf8
set terminal wxt size 900,600 enhanced font "Arial,12"
set grid
set palette rgbformulae 33,13,10
unset colorbox

# Leer encabezado de parámetros
parametros = system("grep '^#' results/poincare.dat | head -n1 | cut -c3-")

# --- Título principal en dos líneas ---
set multiplot layout 1,2 title "Secciones de Poincaré - Sistema Acoplado\n".parametros font "Arial,12"

# ---- Subgráfico 1 ----
set title "(x1, v1)"
set xlabel "x1 (m)"
set ylabel "v1 (m/s)"
plot "results/poincare.dat" using 1:2:1 with points pt 7 ps 0.3 lc palette notitle

# ---- Subgráfico 2 ----
set title "(x2, y2)"
set xlabel "x2 (m)"
set ylabel "v2 (m/s)"
plot "results/poincare.dat" using 3:4:3 with points pt 7 ps 0.3 lc palette notitle

unset multiplot
pause -1 "Presiona Enter para salir"
