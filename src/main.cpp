/**
 * @file     main.cpp
 * @brief    Ocilador de Duffing.
 * @author   Eric Jesús Arciniegas Barreto, Santiago Suarez Sanchez
 * @date     [28/06/2025]
 * @version  1.0
 * @license  Propietario
 *
 */

#include "rungekuta.h"
#include "utilidad.h"

const double alfa = -1;              // término lineal restaurador
const double beta = 3;           // no linealidad cúbica moderada
const std::vector<double> delta = {0.05, 0.05}; // fricción débil (pero no cero)
const double gamma = 1.5;               // fuerza externa moderada
const double omega = 0.6;             // frecuencia de la fuerza externa
const double k = 5;                      // acoplamiento lineal moderado               // frecuencia de la fuerza externa         
const std::vector<double> m = {1.0, 1.0};  // masas iguales
const std::vector<double> l = {0.5, 0.5};   // longitudes





// ——— Declaración de funciones ———

double f1(double, double, double, double, double);
double f2(double, double, double, double, double);

// ——— main ———

int main() {

// Valores iniciales
//double epsilon = 0.01; // pequeña perturbación inicial
std::vector<double> t0 = {0, 70, 0.0001}; // tiempo inicial, final, paso
//std::vector<double> t = {t0.at(0)}, x1 = {M_PI/3}, x2 = {M_PI/4}, y1 = {0}, y2 = {0}; // tiempo, angulos, velociadades angulares  
std::vector<double> t = {t0.at(0)}, x1 = {-1+0.0001}, x2 = {1+0.0001}, y1 = {0}, y2 = {0};

for (int a = 0; t.at(a) < t0.at(1); a++){t.push_back(t0.at(0) + (a + 1) * t0.at(2));}

// Implementación

rk4::rungekuta(t, x1, x2, y1, y2, f1, f2); // llamado de función rungekuta 4
utd::guardado(t, x1, x2, y1, y2,"datos", alfa, beta, gamma, k, omega, m, delta); // llamado de función guardado
rk4::calcularLyapunov("results/datos1.dat","results/datos2.dat","results/lyapunov.dat"); // llamado de función calcularLyapunov
rk4::poincare_grid(-2, 2, 
                   -2, 2,
                   100,
                   1500, 0.001,
                   omega,
                   alfa, beta, gamma,
                   k,
                   m, delta,
                   f1, f2);

utd::animacion1("datos",l,t0.at(2),"anim");
utd::grafica("datos","datos");
}

// ——— Definición de funciones ———

double f1(double t, double x1, double x2, double y1, double y2){
    return -( y1 * delta.at(0) 
            + m.at(0) * alfa * x1                          // término lineal restaurador
            + beta * std::pow(x1,3)                        // no linealidad cúbica Duffing
            + k * ( x1 - x2 )                           // acoplamiento lineal
            + gamma * std::cos(omega * t) )       // fuerza externa periódica
           / m.at(0);
}

double f2(double t, double x1, double x2, double y1, double y2){
    return -( y2 * delta.at(1) 
            + m.at(1) * alfa * x2                          // término lineal restaurador
            + beta * std::pow(x2,3)                        // no linealidad cúbica Duffing
            + k * ( x2 - x1 ) )                         // acoplamiento lineal
           / m.at(1);
}



//double f1(double t, double x1, double x2, double y1, double y2){
//return -( y1 * b.at(0) + m.at(0) * g * std::sin(x1) + k * std::pow(l.at(0),2) * ( x1 - x2 ) + p * std::pow(l.at(0),2) * std::pow(( x1 - x2 ),3) + g0.at(0) * std::pow(l.at(0),2) * std::cos(g0.at(1)*t) ) / ( m.at(0) * std::pow(l.at(0),2) );
//}

//double f2(double t, double x1, double x2, double y1, double y2){
//return -( y2 * b.at(1) + m.at(1) * g * std::sin(x2) + k * std::pow(l.at(1),2) * ( x2 - x1 ) + p * std::pow(l.at(0),2) * std::pow(( x2 - x1 ),3)) / ( m.at(1) * std::pow(l.at(1),2) );
//}

