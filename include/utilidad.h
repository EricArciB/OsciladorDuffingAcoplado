/**
 * @file     utilidad.h
 * @brief    Utilidades para guardado de datos, gráficas dinamicas.
 * @author   Eric Jesús Arciniegas Barreto, Santiago Suarez Sanchez
 * @date     [28/06/2025]
 * @version  1.0
 * @license  Propietario
 *
 * Este paquete contiene las funciones auxiliares para el guardado de datos y gráficados dinamicos (gif).
 */


#ifndef UTILIDAD_H
#define UTILIDAD_H
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>

/// ——— Declaración de funciones ———
/**
 * @brief Defino un nombre de espacios utd por utilidad.
 */
namespace utd{

/**
 * @brief Guarda en un .dat 3 datos guardados en vectores x, y1, y2.
 * @param primer dato a guardar.
 * @param segundo dato a guardar.
 * @param tercer dato a guardar.
 * @param cuarto dato a guardar.
 * @param quinto dato a guardar
 * @param nombre del documento .dat donde se guardaran los datos.
 * @param parámetros del sistema (alfa, beta, gamma, k, omega, m1, m2, delta1, delta2) para incluirlos en el encabezado del archivo de datos.
 */


void guardado(std::vector<double>& t,
              std::vector<double>& x1,
              std::vector<double>& x2,
              std::vector<double>& y1,
              std::vector<double>& y2,
              const std::string& nombreArchivo,
              double alfa, double beta, double gamma,
              double k, double omega,
              const std::vector<double>& m,
              const std::vector<double>& delta);
/**
 * @brief Es una función auxiliar que cuenta las lineas de un documento.
 * @param nombre del documento del que va a contar las lineas.
 */
	int contarLineas(const std::string nombreArchivo);
/**
 * @brief Crea una gráfica dinamica con el uso del script animate.gnu de los pendulos acoplados.
 * @param nombre del documento .dat donde se guardaron los datos.
 * @param longitud de los pendulos.
 * @param double que indica el paso temporal.
 * @param nombre del documento .gif donde se guardaran los datos.
 */
	void animacion1(const std::string nombreArchivo, std::vector<double> l, double dt, const std::string gifSalida);
/**
 * @brief Crea una gráfica de distint.
 * @param nombre del documento .dat donde se guardaron los datos.
 * @param entero de valor 1/0 que indica si se grafica un mapa vectorial predeterminado.
 * @param entero de valor 1/0 que indica si se grafica la trayectoria del .dat.
 * @param nombre del documento .gif donde se guardaran los datos.
 */
	void grafica(const std::string nombreArchivo, const std::string nombreGrafica);
/**
 * @brief Crea una gráfica dinamica con el uso del script animate.gnu, tiene opción de graficar mapas, trayectorias, o ambas en un solo documento de salida.
 * @param nombre del documento .dat donde se guardaron los datos.
 * @param entero de valor 1/0 que indica si se grafica un mapa vectorial predeterminado.
 * @param entero de valor 1/0 que indica si se grafica la trayectoria del .dat.
 * @param nombre del documento .gif donde se guardaran los datos.
 */
	void animacion2(const std::string nombreArchivo, int i,  int j, double dt, const std::string gifSalida);
}

#endif
