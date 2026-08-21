#ifndef DIPOLE_WAVEFUNCTIONS_HPP
#define DIPOLE_WAVEFUNCTIONS_HPP

#include "parameters.hpp"

namespace dipole {

double epsilon2(double z, double Q2, const Parameters& p);

double psiT2(double r, double z, double Q2, const Parameters& p);

double psiL2(double r, double z, double Q2, const Parameters& p);

} // namespace dipole

#endif