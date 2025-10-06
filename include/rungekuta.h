/**
 * @file     rungekuta.h
 * @brief    Método Runge Kuta 4 para 2 odren, con mapeo.
 * @author   Eric Jesús Arciniegas Barreto, Santiago Suarez Sanchez, Esteban Yarik Tobar Diaz.
 * @date     [28/06/2025]
 * @version  1.0
 * @license  Propietario
 *
 * Este paquete contiene la función principal para la implementación del método Runge Kuta 4 para ecuaciones diferenciales de segundo orden, con variables ligadas, incluso en su primera deribada.
 */


#ifndef RUNGEKUTA_H
#define RUNGEKUTA_H
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>



/// ——— Declaración de funciones ———
/**
 * @brief Defino un nombre de espacios Runge-Kuta 4.
 */
namespace rk4{

/**
 * @brief Calcula todos los pasos del método para unos vectores de datos, depende del tamaño del vector parametro t, si t tiene 100 elementos, rungekuta hara 100 pasos, esta preparado para que los pasos no necesariamente sean del mismo tamaño.
 * @param parametro de la variables principales x.
 * @param variable principal que corresponde a la de la primera ecuación diferencial de segundo orden.
 * @param variable principal que corresponde a la de la segunda ecuación diferencial de segundo orden.
 * @param primera variable auxiliar, se uso para caer en dos ecuaciones diferenciales de primer orden.
 * @param segundaa variable auxiliar, se uso para caer en dos ecuaciones diferenciales de primer orden.
 * @param función de la primera variable auxilias, corresponde a su primera derivada.
 * @param función de la segunda variable auxiliar, corresponde a su primera derivada.
 */
	void rungekuta(std::vector<double>& t, std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y1, std::vector<double>& y2,
		       double (*f1)(double, double, double, double, double), double (*f2)(double, double, double, double, double));
/**
 * @brief Calcula un mapa vectorial de la tendencia del ocilador, esto lo hace con una maya del espacio de fase, luego ejecuta un primer paso de rungekuta y usa los valores para calcular el vector de tendencia del punto, imprime un archivo de datos y usa la función gráficar para dibujar el mapa en gnuplot.
 * @param tamaño en x del mapa.
 * @param tamaño en y del mapa.
 * @param el tiempo en el que se quiere graficar el mapa (esto es por si la situacion tiene un espacio de fase dinamico).
 * @param posición en x del mapa.
 * @param posición en y del mapa.
 * @param entero que usara para personalizar el nombre de salida, ejemplot: MapDatos1.dat donde el entero es 1.
 * @param función de la variable principal, corresponde a su primera derivada.
 * @param función de la variable auxiliar, corresponde a su primera derivada.
 */
	void map(double x0, double y0, double x20, double y20, double t0, double a, double b, int i,
		 double (*f1)(double, double, double, double, double), double (*f2)(double, double, double, double, double));
/**
 * @brief Crea un script que gráfica un mapa vectorial en gnuplot, normalizado y ajustado al tamaño del mapa y con esquema de colores para la intencidad de los vectores.
 * @param rango en x.
 * @param rango en y.
 * @param posición en x del mapa (para reajustar el rango).
 * @param posición en y del mapa (para reajustar el rango).
 * @param entero que usara para personalizar el nombre de salida, ejemplot: MapDatos1.gp donde el entero es 1.
 */
	void graficar(double x0, double y0, double a, double b, int i);

void poincare_grid(double y1min, double y1max,
                   double y2min, double y2max,
                   int N,                 // número de divisiones en cada eje
                   double tmax, double h, // tiempo total y paso
                   double omega,          // frecuencia
                   double alfa, double beta, double gamma,
                   double k,
                   const std::vector<double>& m,
                   const std::vector<double>& delta,
                   double (*f1)(double,double,double,double,double),
                   double (*f2)(double,double,double,double,double));
void calcularLyapunov(const std::string &fileA,
                      const std::string &fileB,
                      const std::string &outfile);
}





#endif
