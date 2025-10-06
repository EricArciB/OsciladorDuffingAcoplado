/**
 * @file     map.cpp
 * @brief    Programa auxiliar.
 * @author   Eric Jesús Arciniegas Barreto, Santiago Suarez Sanchez
 * @date     [28/06/2025]
 * @version  1.0
 * @license  Propietario
 *
 * Este script es la base del programa ./bin/auxiliar que se usa para la gráfica dinamica, recibe una entrada que corresponde al tiempo, y con ello actualiza los datos del mapa vectorial predeterminado
 */
 

#include "rungekuta.h"

const double g = 9.81;
const double k = 150;
const double p = 50;
const std::vector<double> g0 = {0,0.5}; // fuerza
const std::vector<double> m = {1,1}; // masa del cuerpo 1 y 2
const std::vector<double> b = {0.5,0.5};// coeficiente de fricción cuerpo 1 y 2
const std::vector<double> l = {0.5,0.5}; // longitud de pendulo del cuerpo 1 y 2

/// ——— Declaración de funciones ———

double f1(double, double, double, double, double);
double f2(double, double, double, double, double);

/// ——— main ———

int main() {
double i;
std::cin >> i;
rk4::map(2, 2.5, 0, 0, i, 0, 0, 0, f1, f2); /// llamar a la función map
std::system("printf \"\\033[A\\033[2K\"");
std::cout << "Tiempo:" << i << " s" << std::endl;
std::system("cp results/MapDatos0.dat scripts/Map.dat");
std::system("rm -f results/MapDatos0.dat maps/MapDatos0.png");
}

/// ——— Definición de funciones ———

double f1(double t, double x1, double x2, double y1, double y2){
return -( y1 * b.at(0) + m.at(0) * g * std::sin(x1) + k * std::pow(l.at(0),2) * ( x1 - x2 ) + p * std::pow(l.at(0),2) * std::pow(( x1 - x2 ),3) + g0.at(0) * std::pow(l.at(0),2) * std::cos(g0.at(1)*t) ) / ( m.at(0) * std::pow(l.at(0),2) );
}

double f2(double t, double x1, double x2, double y1, double y2){
return -( y2 * b.at(1) + m.at(1) * g * std::sin(x2) + k * std::pow(l.at(1),2) * ( x2 - x1 ) + p * std::pow(l.at(0),2) * std::pow(( x2 - x1 ),3)) / ( m.at(1) * std::pow(l.at(1),2) );
}

