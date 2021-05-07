/**
 * @file test_polynome.cpp
 * @brief Programme principal de test pour la mise en oeuvre de la classe polynome
 * @details C'est une batterie de test permettant de s'assure qu'on n'a pas de bogue dans les methode de
 *          la classe polynome
 * @return Toujours le succes !
 */
#include <iostream>
#include "polynome.hpp"

int main()
{
    polynome p({1.,2.,1.});
    std::cout << "p(x) = " << p << std::endl;
    polynome q(3);
    q[0] = 0.; q[1] = 1.; q[2] = 0.; q[3] = -1.;
    std::cout << "q(x) = " << q << std::endl;
    q -= p;
    std::cout << "q(x) -= p(x) => " << q << std::endl;
    polynome r = p - q;
    std::cout << "r(x) = p(x)-q(x) = " << r << std::endl;
    std::cout << "p(x) = " << p << std::endl;
    std::cout << "p(2) = " << p(2) << " ?== 9" << std::endl;
    std::cout << "2.p(x) = " << 2.*p << std::endl;
    return EXIT_SUCCESS;
}
