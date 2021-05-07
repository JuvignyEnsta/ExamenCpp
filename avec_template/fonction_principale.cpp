/**
 * @file fonction_principale.cpp
 * 
 * Programme principal calculant la base orthonormale de polynomes.
 */
#include <vector>
#include <iostream>
#include <cmath>

#include "polynome.hpp"
#include "produits_scalaires.hpp"

/**
 * @brief Construit une base orthonormale de polynôme en utilisant un produit-scalaire fourni par
 *        produit_scalaire_base.
 * @details produit_scalaire_base est une classe abstraite qui permet de définir un opérateur () définissant
 *          produit-scalaire particulier (intégrale avec poids ou non).
 * 
 * @param dim Le nombre de vecteurs orthonormaux à construire 
 * @param dot Le produit scalaire choisi
 * @tparam K  Le type d'éléments utilisé (pour le corps : P appartient à K[X])
 * @return Retourne un ensemble de polynômes orthonormaux pour le produit scalaire passé en argument.
 */
template<typename K>
std::vector<polynome<K>> construire_base_ortho_polynomiale(int dim, const produit_scalaire_base<K>& dot )
{
    std::vector<polynome<K>> base;
    base.reserve(dim);
    // Construit le polynome de degré zéro p0(x) = 1
    base.emplace_back(polynome<K>::canonique(0));
    // Normalisation premier vecteur
    base.back() *= std::sqrt(K(1)/dot(base.back(),base.back()));
    // Algorithme de Gram-Schmidt modifié :
    for ( int i = 1; i < dim; ++i )
    {
        // Construit le polynome de degré i pi(x) = x^i
        polynome<K> q = polynome<K>::canonique(i);
        for ( int j = 0; j < i; ++j )
            q -= dot(q,base[j])*base[j];
        q *= std::sqrt(K(1)/dot(q,q));
        base.emplace_back(q);
    }
    return base;
}

int main()
{
    const double pi = 3.1415926535897932;
    auto dotL = produit_legendre<double>(5);// Produit scalaire int_{-1}^{1}P(x)Q(x)dx
    auto legendre   = construire_base_ortho_polynomiale(5, dotL);
    std::cout << "Legendre : \n";
    for (int i = 0; i < legendre.size(); ++i )
    {
        std::cout << "L_" << i << "(x) = " << legendre[i] << std::endl;
    }
    std::cout << "VERIFICATION\n"
              << "------------\n";
    std::cout << "+1./sqrt(3) doit etre racine de p_2 : ";
    std::cout << "L_2(1./sqrt(3)) = " << legendre[2](1./std::sqrt(3.)) << std::endl;
    std::cout << "-1./sqrt(3) doit etre racine de p_2 : ";
    std::cout << "L_2(-1./sqrt(3)) = " << legendre[2](-1./std::sqrt(3.)) << std::endl;
    std::cout << "+sqrt(3./5.) doit etre racine de p_3 : ";
    std::cout << "L_3(sqrt(3./5.)) = " << legendre[3](std::sqrt(3./5.)) << std::endl;
    std::cout << "-sqrt(3./5.) doit etre racine de p_3 : ";
    std::cout << "L_3(-sqrt(3./5.)) = " << legendre[3](-std::sqrt(3./5.)) << std::endl;
    std::cout << "0 doit etre racine de p_3 : ";
    std::cout << "L_3(0) = " << legendre[3](0) << std::endl;


    auto dot = produit_scalaire_tchebychev<double>(5);// Produit scalaire int_{-1}^{1}P(x)Q(x)/sqrt(1-x²) dx
    auto tchebychev = construire_base_ortho_polynomiale(5, dot);
    std::cout << "Tchebychev : \n";
    std::cout << "T_" << 0 << "(x) = " << tchebychev[0] << std::endl;
    for ( int i = 1; i < tchebychev.size(); ++i )
    {
        std::cout << "T_" << i << "(x) = " << ((1<<(i-1))/tchebychev[i][i])*tchebychev[i] << std::endl;
    }
    std::cout << "VERIFICATION\n"
              << "------------\n";
    std::cout << "Les racines des differents polynomes (a 10^-3 pres a cause de l'integration pas adaptee) : \n";
    for ( int n = 1; n < 5; ++n )
    {
        for ( int i = 1; i <= n; ++i )
        {
            double x = std::cos((2.*i-1.)*pi/(2.*n));
            std::cout << "T_" << n << "(" << x << ") = " << tchebychev[n](x) << std::endl;
        }
    }
    return EXIT_SUCCESS;
}
