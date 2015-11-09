#include <iostream>
#include <iomanip>

int main() {
    double K, G, nu;
    std::cout << "Input the shear modulus (mu1 + mu2):" << std::endl;
    std::cin >> G;
    std::cout << "Input the bulk modulus (K):" << std::endl;
    std::cin >> K;
    nu = K/G;
    nu = (3*nu - 2)/(2*(3*nu + 1));
    std::cout << "Poisson's ratio: " << std::fixed << std::setprecision(8) << nu << std::endl;
}
