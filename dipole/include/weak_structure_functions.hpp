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
// opcional no sentido de que sem LHAPDF o projeto compila, roda e passa
// a bateria de testes -- so nao aceita --use-F3 1.
//
// MAS NAO E DESPREZIVEL EM BAIXA ENERGIA. xF3 e a funcao de estrutura de
// valencia; em E pequeno o x tipico e grande e a valencia domina. Medido
// com o codigo corrigido (GBW, CT10nlo):
//
//     E [GeV]   1e3    1e4    1e5    1e6    1e9    1e14
//     xF3        41%    26%    11%   3.7%  0.87%  0.99%
//
// Ou seja: sem LHAPDF, sigma abaixo de ~1e5 GeV fica subestimada de forma
// significativa. Acima de 1e6 GeV o efeito e de 1% a 4%.
//
// (O numero de 5.8% que circulou antes foi medido no codigo PRE-campanha,
// onde o erro de quadratura inflava F2 por ~100x e fazia o F3 parecer
// irrelevante por comparacao.)
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

    // xF3 de corrente NEUTRA, alvo isoescalar, ordem dominante:
    //
    //     xF3^NC = 2 Soma_q g_Vq g_Aq (x q - x qbar)
    //
    // usando g_L^2 - g_R^2 = g_V g_A. E puramente de VALENCIA, entao cai
    // muito mais rapido com a energia que o termo NC total: acima de
    // ~1e6 GeV contribui menos de 1%. Abaixo de 1e5 GeV nao e
    // desprezivel, pelo mesmo motivo que no CC.
    double xF3_NC_isoscalar(double x, double Q2, BeamType beam) const;

    // true se o binario foi compilado com LHAPDF.
    static bool available();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

BeamType parseBeamType(const std::string& name);

}

#endif
