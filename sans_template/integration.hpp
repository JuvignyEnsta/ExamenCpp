/**
 * @file integration.hpp
 * 
 * Schéma d'intégration pour l'examen (pas forcément adapté à certaines intégrales)
 */
#ifndef _INTEGRATION_HPP_
#define _INTEGRATION_HPP_
#include <string>
#include <cmath>
#include <stdexcept>

using namespace std::string_literals;

/**
 * @brief Formules de quadrature de Gauss-Legendre (pas forcément idéal pour certains produits-scalaires)...
 * @details Calcul l'intégrale entre a et b de la fonction func avec un degré de précision deg à
 *          l'aide de la formule de quadrature de Gauss-Legendre
 * 
 * @param a La borne inférieure de l'intégrale
 * @param b La borne supérieure de l'intégrale
 * @param func La fonction à intégrer
 * @param deg L'ordre d'intégration de la quadrature
 */
template<typename F>
double
quadrature(const double& a, const double& b, const F& func, int deg)
{
    double alpha = 0.5*(b-a);
    double beta  = 0.5*(a+b);
    switch(deg)
    {
    case 1:
        return 2.*alpha*func(beta);
        break;
    case 2:
        {
            double x0 = 1./std::sqrt(3.);
            return alpha*(func(-alpha*x0+beta)+func(alpha*x0+beta));
        }
        break;
    case 3:
        {  
            double x0 = std::sqrt(3./5.);
            double w0 = 5./9.;
            double w1 = 8./9.;
            return alpha*(w0*func(-alpha*x0+beta)+w1*func(beta)+w0*func(alpha*x0+beta));
        }
        break;
    case 4:
        {
            double x0 = std::sqrt(3./7.-2./7.*std::sqrt(6./5.) );
            double w0 = (18.+std::sqrt(30))/36.;
            double x1 = std::sqrt(3./7.+2./7.*std::sqrt(6./5.) );
            double w1 = (18.-std::sqrt(30.) )/36.;
            return alpha*(w0*(func(-alpha*x0+beta)+func(alpha*x0+beta)) + w1*(func(-alpha*x1+beta)+func(alpha*x1+beta)) );
        }
        break;
    case 5:
        {
            double w0 = 128./225.;
            double x1 = std::sqrt(5.-2.*std::sqrt(10./7.))/3.;
            double w1 = (322.+13.*std::sqrt(70.))/900.;
            double x2 = std::sqrt(5.+2.*std::sqrt(10./7.))/3.;
            double w2 = (322.-13.*std::sqrt(70.))/900.;
            return alpha*(w0*func(beta) + w1*(func(-alpha*x1+beta)+func(alpha*x1+beta) ) +
                          w2*(func(-alpha*x2+beta)+func(alpha*x2+beta)) );
        }
        break;
    default:
        std::string error = "Pas de quadrature pour l'ordre "s + std::to_string(deg) + "."s;
        throw std::runtime_error(error);
    }
}

/**
 * @brief Calcul de l'intégrale de la fonction func entre les bornes a et b
 * @details Calcul de l'intégrale de la fonction func entre les bornes a et b en subdivisant
 *          l'intervalle [a,b] en sous-intervalles
 * 
 * @param a La borne inférieure de l'intégrale
 * @param b La borne supérieure de l'intégrale
 * @param func La fonction à intégrer
 * @param deg L'ordre de quadrature pour intégrer sur chaque sous-intervalle.
 */
template<typename F>
double
integration(const double& a, const double& b, const F& func, int deg)
{
    const int N = 100000;
    double val = 0.;
    double h   = (b-a)/N;
    for ( int i = 0; i < N; ++i )
    {
        val += quadrature(a+i*h, a + (i+1)*h, func, deg);
    }
    return val;
}
#endif
