#include <iostream>
#include <cmath>

double thomas_fermi_screening(double n, double Z, double r) {
    double a0 = 0.529; // Bohr radius in Angstroms
    double e = 1.602e-19; // Elementary charge in Coulombs
    double hbar = 1.055e-34; // Reduced Planck constant in Joule-seconds
    double me = 9.109e-31; // Electron mass in kilograms
    double eps0 = 8.854e-12; // Vacuum permittivity in Farads per meter
    double pi = 3.14159;

    double rs = pow(3.0/(4.0*pi*n), 1.0/3.0); // Wigner-Seitz radius
    double kf = pow(3.0*pi*pi*n, 1.0/3.0)/a0; // Fermi wavevector
    double ef = pow(hbar*kf, 2.0)/(2.0*me); // Fermi energy
    double eps = 1.0 + rs*(Z/a0)*sqrt(ef/(2.0*pi*eps0)); // Dielectric constant
    double vtf = -e*e/(4.0*pi*eps0*r)*exp(-kf*r)/(1.0+kf*r); // Thomas-Fermi potential
    double vext = 0.0; // External potential

    return vtf + vext;
}

int main() {
    double n = 1.0e20; // Electron density in m^-3
    double Z = 29.0; // Atomic number of copper
    double r = 1.0e-10; // Distance from nucleus in meters

    double vtf = thomas_fermi_screening(n, Z, r);
    std::cout << "Thomas-Fermi potential: " << vtf << " J" << std::endl;

    return 0;
}