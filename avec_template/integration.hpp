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
template<typename K, typename F>
K
quadrature(const K& a, const K& b, const F& func, int deg)
{
    double alpha = K(0.5)*(b-a);
    double beta  = K(0.5)*(a+b);
    switch(deg)
    {
    case 1:
        return K(2)*alpha*func(beta);
        break;
    case 2:
        {
            K x0 = K(1)/std::sqrt(K(3));
            return alpha*(func(-alpha*x0+beta)+func(alpha*x0+beta));
        }
        break;
    case 3:
        {  
            K x0 = std::sqrt(K(3)/K(5));
            K w0 = K(5)/K(9);
            K w1 = K(8)/K(9);
            return alpha*(w0*func(-alpha*x0+beta)+w1*func(beta)+w0*func(alpha*x0+beta));
        }
        break;
    case 4:
        {
            K x0 = std::sqrt(K(3)/K(7)-K(2)/K(7)*std::sqrt(K(6)/K(5)) );
            K w0 = (K(18)+std::sqrt(K(30)) )/K(36);
            K x1 = std::sqrt(K(3)/K(7)+K(2)/K(7)*std::sqrt(K(6)/K(5)) );
            K w1 = (K(18)-std::sqrt(K(30)) )/K(36);
            return alpha*(w0*(func(-alpha*x0+beta)+func(alpha*x0+beta)) + w1*(func(-alpha*x1+beta)+func(alpha*x1+beta)) );
        }
        break;
    case 5:
        {
            K w0 = K(128)/K(225);
            K x1 = std::sqrt(K(5)-K(2)*std::sqrt(K(10)/K(7)) )/K(3);
            K w1 = (K(322)+K(13)*std::sqrt(70) )/K(900);
            K x2 = std::sqrt(K(5)+K(2)*std::sqrt(K(10)/K(7)) )/K(3);
            K w2 = (K(322)-K(13)*std::sqrt(70) )/K(900);
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
template<typename K, typename F>
K
integration(const K& a, const K& b, const F& func, int deg)
{
    const int N = 100000;
    K val = 0.;
    K h   = (b-a)/N;
    for ( int i = 0; i < N; ++i )
    {
        val += quadrature(a+i*h, a + (i+1)*h, func, deg);
    }
    return val;
}
#endif
