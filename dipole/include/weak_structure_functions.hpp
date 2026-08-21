#ifndef WEAK_STRUCTURE_FUNCTIONS_HPP
#define WEAK_STRUCTURE_FUNCTIONS_HPP

#include <memory>
#include <string>
#include <LHAPDF/LHAPDF.h>

namespace weak {

enum class BeamType {
    Neutrino,
    AntiNeutrino
};

class WeakStructureFunctions {
public:
    explicit WeakStructureFunctions(
        const std::string& pdf_set = "NNPDF31_nlo_as_0118",
        int member = 0
    );

    double xF3_CC_isoscalar(double x, double Q2, BeamType beam) const;

private:
    std::unique_ptr<LHAPDF::PDF> pdf_;

    double xf(int pid, double x, double Q2) const;
};

BeamType parseBeamType(const std::string& name);

}

#endif