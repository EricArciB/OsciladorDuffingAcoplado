##
## \file animate.gnu
## \brief Script de Gnuplot para generar una animación del espacio de fase.
## \details Este script produce un GIF animado a partir de datos generados por un programa de c++. 
##          Los datos representan la evolución de pendulos acoplados por un resorte.

# Terminal animado en formato GIF con retardo entre cuadros
set terminal gif animate delay 1 size 800,600

# Archivo de salida
set output 'results/animate1.gif'

# Rango de los ejes
set xrange [-1.1*(L1+L2):1.1*(L1+L2)]
set yrange [-1.2*(L1+L2)/2:0.1*(L1+L2)]

# Etiquetas
set xlabel 'Posición x (m)'
set ylabel 'Posición y (m)'

# Cuadrícula
set grid


# Bucle para generar cada cuadro de la animación
do for [i=1:N] {

	# Tiempo:
	set title sprintf('Pendulos Acoplados Tiempo: %.2f s', i*T)
	
	# Cuerdas de pendulo:
	
	cmd1 = sprintf("awk 'NR==%d {print $2}' results/animate1.dat", i)
	cmd2 = sprintf("awk 'NR==%d {print $3}' results/animate1.dat", i)
	
	theta1 = real(system(cmd1))
	theta2 = real(system(cmd2))
	 
	l11 = -(L1+L2)/2
	l12 = 0
	l13 = -(L1+L2)/2 + L1*sin(theta1)
	l14 = -L1*cos(theta1)
	
	set print 'scripts/l1.dat'
	print l11, l12
	print l13, l14
	unset print
	
	l21 = (L1+L2)/2
	l22 = 0
	l23 = (L1+L2)/2 + L2*sin(theta2)
	l24 = -L2*cos(theta2)
	
	set print 'scripts/l2.dat'
	print l21, l22
	print l23, l24
	unset print
	
	# Resorte
	
	# Puntos inicial y final
	x1 = -(L1+L2)/2 + L1*sin(theta1)
	y1 = - L1*cos(theta1)
	x2 = (L1+L2)/2 + L2*sin(theta2)
	y2 = - L2*cos(theta2)

	# Vector dirección y longitud
	dx = x2 - x1
	dy = y2 - y1
	L = sqrt(dx**2 + dy**2)
	
	# Parámetros del seno
	A = 0.2*(L1+L2)/2
	k = 2*pi*15 / L  # dos ciclos

	# Unitarios
	ux = dx / L
	uy = dy / L
	nx = -dy / L
	ny = dx / L
	
	# Funciones para resorte
	set parametric
	set trange [0:L]
	set size ratio -1
	unset key

	x(t) = x1 + t*ux + A*sin(k*t)*nx
	y(t) = y1 + t*uy + A*sin(k*t)*ny

	plot 'scripts/l1.dat' using 1:2 w l lc rgb 'black' notitle, \
	     'scripts/l2.dat' using 1:2 w l lc rgb 'black' notitle, \
	     'results/animate1.dat' every ::i-1::i-1 using (-(L1+L2)/2 + L1*sin($2)):(-L1*cos($2)) w p pt 7 ps 3 title 'm1', \
	     'results/animate1.dat' every ::i-1::i-1 using ((L1+L2)/2 + L2*sin($3)):(-L2*cos($3)) w p pt 7 ps 3 title 'm2', \
	     x(t), y(t) with lines lw 2 lc rgb 'blue'

}

# Termina el archivo de salida
unset output
