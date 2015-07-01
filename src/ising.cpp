#include <stddef.h>
#include <iostream>
#include <gsl/gsl_sf_bessel.h>

#include "Lattice.h"

int main(int argc, char** argv) {
    Lattice* lat = new Lattice(64);
    std::cout << gsl_sf_bessel_J0(0.) << std::endl;
}
