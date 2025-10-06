##
## \file animate.gnu
## \brief Script de Gnuplot para generar una animación del espacio de fase.
## \details Este script produce un GIF animado a partir de datos generados por un programa auxiliar. 
##          Los datos representan la evolución de un sistema dinámico en el espacio de fase (posición vs velocidad).
##          Se grafican vectores normalizados con color proporcional a su módulo.

# Terminal animado en formato GIF con retardo entre cuadros
set terminal gif animate delay 1 size 800,600

# Archivo de salida
set output 'results/animate2.gif'

# Rango de los ejes
set xrange [-2:2]
set yrange [-2.5:2.5]

# Etiquetas y título
set title 'Espacio de fase'
set xlabel 'Posición (m)'
set ylabel 'Velocidad (m/s)'

# Cuadrícula
set grid

# Paleta de colores y rango de la barra de color
set palette rgb 33,13,10
set cbrange [0:14]
set colorbox

# Bucle para generar cada cuadro de la animación
do for [i=1:N] {

	# Ejecuta el programa auxiliar con el tiempo correspondiente
	system sprintf("echo %.2f | ./bin/auxiliar", i*T)

	# Grafica:
	# - Vectores de velocidad normalizados con color según módulo.
	# - Línea que muestra la evolución temporal hasta el tiempo i*T.
	plot 'scripts/Map.dat' using (C ? $1 : 1/0):2:(0.1*($3/sqrt($3**2+$4**2))):(0.1*($4/sqrt($3**2+$4**2))):(sqrt($3**2+$4**2)) with vectors head filled lc palette notitle, 'results/animate2.dat' every ::1::i using (D ? $2 : 1/0):3 with line title sprintf('Tiempo: %.2f s', i*T)
}

# Termina el archivo de salida
unset output
