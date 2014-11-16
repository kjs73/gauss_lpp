#include <iostream>
#include <stdexcept>
#include <cstdlib>

#include "schramm_equation.hpp"

int main(int argc, char* argv[])
{
    if (argc != 3) {
        throw std::runtime_error("wrong number of arguments");
    }
    const double kappa = atof(argv[1]);
    const double phi = atof(argv[2]);
    std::cout.precision(std::numeric_limits<double>::digits10);
    std::cout << schramm_equation(kappa)(phi) << "\n";
    return 0;
}
