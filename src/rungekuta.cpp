/**
 * @file     rungekuta.cpp
 * @brief    Método Runge Kuta 4 para 2 odren, con mapeo.
 * @author   Eric Jesús Arciniegas Barreto, Santiago Suarez Sanchez, Esteban Yarik Tobar Diaz.
 * @date     [28/06/2025]
 * @version  1.0
 * @license  Propietario
 *
 * Este archivo contiene la implementación de las funciones declaradas en rungekuta.h
 */

#include "rungekuta.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>

namespace rk4{

/// -----------------------------------------------------------------------------
/// 1. Método Runge Kuta 4 para ecuación diferencial de segundo orden.
/// -----------------------------------------------------------------------------
	void rungekuta(std::vector<double>& t, std::vector<double>& x1, std::vector<double>& x2, std::vector<double>& y1, std::vector<double>& y2,
		       double (*f1)(double, double, double, double, double), double (*f2)(double, double, double, double, double)){
		int j = 0; double h; std::vector<double> k1(4), k2(4), l1(4), l2(4);
		
		for (std::size_t i=0; i+1 < t.size(); i++){
			h = t.at(i+1) - t.at(i);
			
			k1.at(0) = h * y1.at(i);
			k2.at(0) = h * y2.at(i);
			k1.at(1) = h * ( y1.at(i) + k1.at(0)/2 );
			k2.at(1) = h * ( y2.at(i) + k2.at(0)/2 );
			k1.at(2) = h * ( y1.at(i) + k1.at(1)/2 );
			k2.at(2) = h * ( y2.at(i) + k2.at(1)/2 );
			k1.at(3) = h * ( y1.at(i) + k1.at(2) );
			k2.at(3) = h * ( y2.at(i) + k2.at(2) );
			
			l1.at(0) = h * f1(t.at(i), x1.at(i), x2.at(i), y1.at(i), y2.at(i));
			l2.at(0) = h * f2(t.at(i), x1.at(i), x2.at(i), y1.at(i), y2.at(i));
			l1.at(1) = h * f1(t.at(i) + h/2, x1.at(i) + k1.at(0)/2, x2.at(i) + k2.at(0)/2,
							 y1.at(i) + l1.at(0)/2, y2.at(i) + l2.at(0)/2);
			l2.at(1) = h * f2(t.at(i) + h/2, x1.at(i) + k1.at(0)/2, x2.at(i) + k2.at(0)/2,
							 y1.at(i) + l1.at(0)/2, y2.at(i) + l2.at(0)/2);
			l1.at(2) = h * f1(t.at(i) + h/2, x1.at(i) + k1.at(1)/2, x2.at(i) + k2.at(1)/2,
							 y1.at(i) + l1.at(1)/2, y2.at(i) + l2.at(1)/2);
			l2.at(2) = h * f2(t.at(i) + h/2, x1.at(i) + k1.at(1)/2, x2.at(i) + k2.at(1)/2,
							 y1.at(i) + l1.at(1)/2, y2.at(i) + l2.at(1)/2);
			l1.at(3) = h * f1(t.at(i) + h, x1.at(i) + k1.at(2), x2.at(i) + k2.at(2), y1.at(i) + l1.at(2), y2.at(i) + l2.at(2));
			l2.at(3) = h * f2(t.at(i) + h, x1.at(i) + k1.at(2), x2.at(i) + k2.at(2), y1.at(i) + l1.at(2), y2.at(i) + l2.at(2));
			
			x1.push_back(x1.at(i) + (k1.at(0) + 2*k1.at(1) + 2*k1.at(2) + k1.at(3))/6);
			x2.push_back(x2.at(i) + (k2.at(0) + 2*k2.at(1) + 2*k2.at(2) + k2.at(3))/6);
			y1.push_back(y1.at(i) + (l1.at(0) + 2*l1.at(1) + 2*l1.at(2) + l1.at(3))/6);
			y2.push_back(y2.at(i) + (l2.at(0) + 2*l2.at(1) + 2*l2.at(2) + l2.at(3))/6);	


			if (!std::isfinite(x1.back()) || !std::isfinite(y1.back()) ||
    !std::isfinite(x2.back()) || !std::isfinite(y2.back())) {
    std::cerr << "NaN detectado en paso i=" << i
              << " t=" << t.at(i)
              << " x1=" << x1.back() << " y1=" << y1.back()
              << " x2=" << x2.back() << " y2=" << y2.back() << "\n";
    return; // o break
}

			j += 1;	
		}
		
		//    std::cout<<"Iteraciones: "<< j <<"\n";
	}
/// -----------------------------------------------------------------------------
/// 2. Mapeo vectorial del espacio de fase.
/// -----------------------------------------------------------------------------
	void map(double x0, double y0, double x20, double y20, double t0, double a, double b, int i,
	 double (*f1)(double, double, double, double, double), double (*f2)(double, double, double, double, double)){
		std::vector<double> t = {t0, t0 + 0.01}, x1, y1, x2 = {x20}, y2 = {y20}, v;
		std::ofstream data("results/MapDatos" + std::to_string(i) + ".dat");
		if(data.is_open()){
		data.precision(5);
		data << std::fixed;
			for(int j = 1; j < 21; j++){
			for(int k = 1; k < 21; k++){
			x1.push_back(a + std::abs(1.05*x0) - j * std::abs(x0/10)); y1.push_back(b + std::abs(1.05*y0) - k * std::abs(y0/10));
			rungekuta(t,x1,x2,y1,y2,f1,f2);
			double m = (y1.at(1)-y1.at(0))/(x1.at(1)-x1.at(0)); /// Pendiente de tendencia, se usara para calcular el vector
			v.push_back(std::abs(m) * std::cos(std::atan(m)));
			v.push_back(std::abs(m) * std::sin(std::atan(m)));
			data << x1.at(0) << "\t" << y1.at(0) << "\t" << v.at(0) << "\t" << v.at(1) << "\n";
			x1.clear();y1.clear();v.clear();
			}
			}
		data.close();
		graficar(x0,y0,a,b,i);
		std::string comando = "gnuplot -p scripts/MapDatos" + std::to_string(i) + ".gp";
		std::system(comando.c_str()); /// Ejecuta el script para gráficar el mapa
		comando = "rm -f scripts/MapDatos" + std::to_string(i) + ".gp";
		std::system(comando.c_str()); /// Elimina el archivo .gp 
		}else{
		std::cerr << "Error: No se pudo abrir el archivo para escritura.\n";}
	}
/// -----------------------------------------------------------------------------
/// 3. Generador de gráfica para mapeo vectorial.
/// -----------------------------------------------------------------------------
	void graficar(double x0, double y0, double a, double b, int i){
		std::ofstream script("scripts/MapDatos" + std::to_string(i) + ".gp");
		script.precision(5);
		script << std::fixed;
		if(script.is_open()){
		script << "set xrange [" << a - std::abs(x0) << ":" << a + std::abs(x0) <<"]\n";
		script << "set yrange [" << b - std::abs(y0) << ":" << b + std::abs(y0) <<"]\n";	
		script << "set grid\n";
		script << "set palette rgb 33,13,10\n";
		script << "set cbrange [0:14]\n";
		script << "set colorbox\n";
		script << "set term pngcairo size 800,600 enhanced\n";
		script << "set output 'maps/MapDatos"<< i <<".png'\n";
		script << "plot 'results/MapDatos"<< i <<".dat' using 1:2:(";
		script << "0.5 *" << std::abs(x0) << "/10 *($3/sqrt($3**2+$4**2))):(";
		script << "0.5 *" << std::abs(x0) << "/10 *($4/sqrt($3**2+$4**2))):(";
		script << "sqrt($3**2+$4**2)) with vectors head filled lc palette notitle\n";
		script << "unset output\n";
		}else{
		std::cerr << "Error: No se pudo abrir el script para escritura.\n";}
	}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Poincaré grid para Duffing acoplados
// Genera una grilla de condiciones iniciales y guarda (x1,y1,x2,y2)
// en t = n*T, con T = 2*pi/omega

void poincare_grid(double y1min, double y1max,
                   double y2min, double y2max,
                   int N,
                   double tmax, double h,
                   double omega,
                   double alfa, double beta, double gamma,
                   double k,
                   const std::vector<double>& m,
                   const std::vector<double>& delta,
                   double (*f1)(double,double,double,double,double),
                   double (*f2)(double,double,double,double,double))
{
    double T = 2.0 * 3.14159265358979323846 / omega;
    int nsteps = static_cast<int>(tmax / h);

    std::ofstream fout("results/poincare.dat");
    if (!fout.is_open()) {
        std::cerr << "Error al abrir results/poincare.dat\n";
        return;
    }

    fout << std::fixed << std::setprecision(6);

    // Encabezado con parámetros del sistema 
fout.unsetf(std::ios::fixed); // desactiva formato fijo
fout << std::setprecision(6); // deja precisión general (solo muestra lo necesario)
fout << "# alfa=" << alfa
     << ", beta=" << beta
     << ", gamma=" << gamma
     << ", k=" << k
     << ", omega=" << omega
     << ", m1=" << m[0]
     << ", m2=" << m[1]
     << ", delta1=" << delta[0]
     << ", delta2=" << delta[1]
     << std::endl;
fout << std::fixed << std::setprecision(6); // vuelve al formato fijo para los datos numéricos


    // recorrer grilla
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {

            double x1_0 = y1min + i * (y1max - y1min) / (N - 1); 
            double x2_0 = y2min + j * (y2max - y2min) / (N - 1);
            double y1_0 = 1;//-3.0 + 6.0 * i / (N-1); // rango [-2,2]
            double y2_0 = 1;//-3.0 + 6.0 * j / (N-1); // rango [-2,2]


            std::vector<double> t, x1, y1, x2, y2;
            t.push_back(0.0);
            x1.push_back(x1_0);
            y1.push_back(y1_0);
            x2.push_back(x2_0);
            y2.push_back(y2_0);

            for (int k = 1; k <= nsteps; k++) t.push_back(k * h);

            rungekuta(t, x1, x2, y1, y2, f1, f2);

            int Nskip = 35;
            for (int n = Nskip; n * T <= tmax; n++) {
                double t_target = n * T;
                size_t k = static_cast<size_t>(t_target / h);
                if (k < t.size()) {
                    fout << x1[k] << "\t" << y1[k] << "\t"
                         << x2[k] << "\t" << y2[k] << "\n";
                }
            }
        }
    }

    fout.close();
}

// -----------------------------------------------------------------------------//
// ------------------------------------------------------------------------------
// Cálculo del exponente de Lyapunov a partir de dos archivos de datos
// -----------------------------------------------------------------------------

struct Punto {
    double t, x1, x2, y1, y2;
};

void calcularLyapunov(const std::string &fileA,
                      const std::string &fileB,
                      const std::string &outfile) {
    std::ifstream file0(fileA);
    std::ifstream file1(fileB);

    if (!file0.is_open() || !file1.is_open()) {
        std::cerr << "Error abriendo archivos: "
                  << fileA << " o " << fileB << std::endl;
        return;
    }

    // Saltar encabezado
    std::string header;
    getline(file0, header);
    getline(file1, header);

    std::vector<Punto> traer0, traer1;
    Punto p0, p1;

    // Leer datos de ambos archivos
    while (file0 >> p0.t >> p0.x1 >> p0.x2 >> p0.y1 >> p0.y2 &&
           file1 >> p1.t >> p1.x1 >> p1.x2 >> p1.y1 >> p1.y2) {
        traer0.push_back(p0);
        traer1.push_back(p1);
    }

    file0.close();
    file1.close();

    int N = traer0.size();
    if (N < 2) {
        std::cerr << "Muy pocos datos para calcular Lyapunov" << std::endl;
        return;
    }

    // Distancia inicial
    double dx1 = traer0[0].x1 - traer1[0].x1;
    double dx2 = traer0[0].x2 - traer1[0].x2;
    double dy1 = traer0[0].y1 - traer1[0].y1;
    double dy2 = traer0[0].y2 - traer1[0].y2;
    double d0 = sqrt(dx1*dx1 + dx2*dx2 + dy1*dy1 + dy2*dy2);

    if (d0 <= 1e-12) d0 = 1e-12; // evita log(0)

    std::ofstream fout(outfile);
    if (!fout.is_open()) {
        std::cerr << "No se pudo crear " << outfile << std::endl;
        return;
    }

    fout << std::fixed << std::setprecision(6);
    fout << std::setw(10) << "t(s)"
         << std::setw(20) << "distancia"
         << std::setw(20) << "ln(d/d0)"
         << std::setw(20) << "lambda(t)" << std::endl;

    double suma_lambda = 0.0;
    int contador = 0;

    for (int i = 1; i < N; i++) {
        double dx1 = traer0[i].x1 - traer1[i].x1;
        double dx2 = traer0[i].x2 - traer1[i].x2;
        double dy1 = traer0[i].y1 - traer1[i].y1;
        double dy2 = traer0[i].y2 - traer1[i].y2;

        double dist = sqrt(dx1*dx1 + dx2*dx2 + dy1*dy1 + dy2*dy2);
        if (dist <= 1e-12 || traer0[i].t <= 1e-12) continue; // evita log(0) y división por 0

        double ln_ratio = log(dist / d0);
        double lambda_t = ln_ratio / traer0[i].t;

        suma_lambda += lambda_t;
        contador++;

        fout << std::setw(10) << traer0[i].t
             << std::setw(20) << dist
             << std::setw(20) << ln_ratio
             << std::setw(20) << lambda_t << std::endl;
    }

    fout.close();

    if (contador > 0) {
        double lambda_prom = suma_lambda / contador;
        std::cout << "Lambda promedio = " << lambda_prom << std::endl;
    } else {
        std::cout << "No se pudo calcular lambda promedio." << std::endl;
    }

    std::cout << "Lyapunov guardado en " << outfile << std::endl;
}

// -----------------------------------------------------------------------------




}


