#ifndef WEAK_STRUCTURE_FUNCTIONS_HPP
#define WEAK_STRUCTURE_FUNCTIONS_HPP

#include <memory>
#include <string>

namespace weak {

enum class BeamType {
    Neutrino,
    AntiNeutrino
};

// xF3 de corrente carregada, a partir de PDFs colineares do LHAPDF.
//
// E a UNICA parte do projeto que precisa de biblioteca externa, e e
// opcional: sem LHAPDF o projeto compila e roda normalmente, so nao
// aceita --use-F3 1. A contribuicao de xF3 e de valencia, entao decai
// depressa com a energia (5.8% em 1e3 GeV, 0.006% em 1e9).
//
// O tipo do LHAPDF fica escondido num pimpl para que este cabecalho
// nao dependa dele.
class WeakStructureFunctions {
public:
    explicit WeakStructureFunctions(
        const std::string& pdf_set = "CT10nlo",
        int member = 0
    );
    ~WeakStructureFunctions();

    WeakStructureFunctions(const WeakStructureFunctions&) = delete;
    WeakStructureFunctions& operator=(const WeakStructureFunctions&) = delete;

    double xF3_CC_isoscalar(double x, double Q2, BeamType beam) const;

    // true se o binario foi compilado com LHAPDF.
    static bool available();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

BeamType parseBeamType(const std::string& name);

}

#endif
