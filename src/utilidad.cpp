/**
 * @file     utilidad.cpp
 * @brief    Utilidades para guardado de datos, gráficas dinamicas.
 * @author   Eric Jesús Arciniegas Barreto, Santiago Suarez Sanchez
 * @date     [28/06/2025]
 * @version  1.0
 * @license  Propietario
 *
 * Este archivo contiene la implementación de las funciones declaradas en utilidad.h
 */


#include "utilidad.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
namespace utd{
/// -----------------------------------------------------------------------------
/// 1. Guardado de datos.
/// -----------------------------------------------------------------------------


void guardado(std::vector<double>& t,
              std::vector<double>& x1,
              std::vector<double>& x2,
              std::vector<double>& y1,
              std::vector<double>& y2,
              const std::string& nombreArchivo,
              double alfa, double beta, double gamma,
              double k, double omega,
              const std::vector<double>& m,
              const std::vector<double>& delta)
{
    std::ofstream fout("results/" + nombreArchivo + ".dat");

    if (!fout.is_open()) {
        std::cerr << "Error al abrir el archivo results/" << nombreArchivo << ".dat\n";
        return;
    }

    // Encabezado con parámetros del sistema 
    fout.unsetf(std::ios::fixed);
    fout << std::setprecision(6);
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

    // Datos simulados
    for (std::size_t i = 0; i < t.size(); ++i) {
        fout << t.at(i) << "\t"<< x1.at(i) << "\t"<< x2.at(i) << "\t"<< y1.at(i) << "\t"<< y2.at(i) << "\n";
    }

    fout.close();
}


/// -----------------------------------------------------------------------------
/// 2. Contado de lineas de un archivo.
/// -----------------------------------------------------------------------------
	int contarLineas(const std::string nombreArchivo){
	
		std::ifstream archivo(nombreArchivo);
		
		int lineas = 0; std::string linea;
		while (getline(archivo, linea)){lineas++;}
		
    		archivo.close();
    	return lineas;
	}
/// -----------------------------------------------------------------------------
/// 3. Generar una animación.
/// -----------------------------------------------------------------------------
	void animacion1(const std::string nombreArchivo, std::vector<double> l, double dt, const std::string gifSalida){
	
		std::string instru1 = "cp results/" + nombreArchivo + ".dat results/animate1.dat"; 	
		std::string instru2 = "gnuplot -e \"L1=" + std::to_string(l.at(0)) + "; L2=" + std::to_string(l.at(1)) + "; T="+ std::to_string(dt) + "; N=" + std::to_string(contarLineas("results/"+ nombreArchivo +".dat")) + "\" -p scripts/animate1.gnu";
		std::string instru3 = "cp results/animate1.gif results/" + gifSalida + ".gif";
		std::system(instru1.c_str());
		std::system(instru2.c_str());
		std::system(instru3.c_str());
		std::system("rm -f results/animate1.gif results/animate1.dat");
	}
/// -----------------------------------------------------------------------------
/// 4. Generar graficas.
/// -----------------------------------------------------------------------------
	void grafica(const std::string nombreArchivo, const std::string nombreGrafica){
		std::string instru1 = "cp results/" + nombreArchivo + ".dat results/grafica.dat"; 	
		std::string instru2 = "gnuplot -p scripts/graficas.gp";
		std::string instru3 = "cp results/grafica.eps results/" + nombreGrafica + ".eps";
		std::system(instru1.c_str());
		std::system(instru2.c_str());
		std::system(instru3.c_str());
		std::system("rm -f results/grafica.eps results/grafica.dat");
	
	}
/// -----------------------------------------------------------------------------
/// 5. Generar una animación.
/// -----------------------------------------------------------------------------
	void animacion2(const std::string nombreArchivo, int i,  int j, double dt, const std::string gifSalida){
		std::string instru1 = "cp results/" + nombreArchivo + ".dat results/animate2.dat"; 	
		std::string instru2 = "gnuplot -e \"C=" + std::to_string(i) + "; D=" + std::to_string(j) + "; T="+ std::to_string(dt) + "; N=" + std::to_string(contarLineas("results/"+ nombreArchivo +".dat")) + "\" -p scripts/animate2.gnu";
		std::string instru3 = "cp results/animate2.gif results/" + gifSalida + ".gif";
		std::system(instru1.c_str());
		std::system(instru2.c_str());
		std::system(instru3.c_str());
		std::system("rm -f results/animate2.gif results/animate2.dat");
	}
}


