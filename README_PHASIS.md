# PHASIS — transporte DIS ao longo de geodésicas nulas

Profundidade óptica de neutrinos de ultra-alta energia ao longo de geodésicas
nulas em espaço-tempo estático e esfericamente simétrico, com perfil de matéria
arbitrário.

**Fase 1 (esta):** espaço plano. Sem relatividade geral, sem disco, sem
regeneração. O objetivo é uma base validada contra resultados analíticos antes
de complicar.

## Unidades

Regra única, sem exceção — dentro do código tudo está em:

| grandeza | unidade |
|---|---|
| comprimento | cm |
| densidade | g/cm³ |
| seção de choque | cm² |
| energia | GeV |
| ângulo | rad |

Constantes de conversão (`include/phasis/units.hpp`) existem para **entrada e
saída**, nunca para uso interno. A checagem dimensional de τ deve fechar a olho:

```
tau = [mol⁻¹] · [g/cm³] · [cm²] · [cm] · [mol/g]   →   adimensional
```

## Física

Métrica `ds² = -f(r)dt² + h(r)dr² + r²dΩ²`. Geodésica nula com parâmetro de
impacto `b = L/E_inf`.

Elemento de comprimento próprio medido pelo observador estático:

```
dl/dr = √h(r) / √(1 - f(r) b²/r²)
```

Energia local: `E_loc(r) = E_inf/√f(r)`. Em Minkowski é a identidade, mas já
está escrita assim — na Fase 2 o redshift faz σ variar ao longo do raio, e
nada aqui precisa mudar. **σ fica dentro da integral**, mesmo quando fatoraria.

```
tau = ∫ N_A · ρ(r,θ) · σ(E_loc(r)) · dl        P = exp(-tau)
```

## Duas armadilhas numéricas, e como são tratadas

**1. A singularidade no ponto de retorno.** O integrando vai como
`(r - r_t)^(-1/2)`. A substituição `r = r_t + s²` (`dr = 2s ds`) a cancela
exatamente.

**2. O cancelamento catastrófico, que é menos óbvio.** Calcular
`w = 1 - f b²/r²` diretamente perto de `r_t` subtrai dois números quase iguais:
em `r - r_t ~ 10⁻¹⁶ r` não sobra dígito. Usando `b² = r_t²/f(r_t)`:

```
w(r) = (r - r_t)·T(r) / (f(r_t) r²),    T(r) = [f(r_t)r² - f(r)r_t²]/(r - r_t)
```

`T` é a parte suave. Com `r - r_t = s²` o fator `s` **cancela analiticamente** e
o integrando vira

```
dl/dr · dr/ds = 2r·√( h(r)·f(r_t)/T(r) )
```

sem singularidade e sem cancelamento. Cada métrica fornece `T` em forma fechada
(`Metric::turning_factor`); em Minkowski, `T = r + r_t`, exato. É isso que dá os
10⁻¹⁶ dos testes T1/T2/T4 — não a tolerância do integrador.

**Os dois ramos** (entrada e saída) são integrados separadamente. Em perfil
esférico eles são iguais, mas o código não assume isso: na Fase 3 `ρ` depende de
`θ` e a simetria quebra.

## Arquitetura

`trace_ray` é uma **função pura** — sem estado global, sem I/O, nada estático.
Um `#pragma omp parallel for` por cima dela é legítimo sem nenhuma mudança.

```cpp
Result trace_ray(const Ray&, const Metric&, const DensityProfile&,
                 const CrossSection&, const IntegratorOpts&);
```

Implementações desta fase: `Minkowski`, `UniformBall`, `PowerLawHalo`,
`PowerLawCrossSection`, `TableCrossSection`.

## Consumindo a produção do `dipole`

O acoplamento entre os dois subprojetos é **por arquivo**, não por link:
`TableCrossSection` lê tanto CSV `E_GeV,sigma_cm2` quanto a saída de três
colunas do `dipole` (`Enu_GeV  sigma_GeV_minus2  sigma_cm2`), e preserva o
cabeçalho `#` de proveniência.

**Não extrapola.** Fora da faixa tabelada lança `std::out_of_range` com a faixa
no texto do erro. Estender uma seção de choque de saturação para fora do
intervalo em que foi calculada é exatamente o tipo de coisa que produz resultado
errado sem aviso.

## Compilar e rodar

```bash
make            # ou: mkdir build && cd build && cmake .. && make
make test       # suite de aceitação T1..T6

./build/trace data/exemplo_terra.cfg
```

`./build/trace` sem argumento imprime as chaves de configuração aceitas.

## Testes de aceitação — Fase 1

| | o que verifica | resultado |
|---|---|---|
| T1 | esfera homogênea, b=0: `τ = N_A ρ₀ σ 2R` | 3,7×10⁻¹⁶ |
| T2 | b/R ∈ {0; 0,1; 0,5; 0,9; 0,99}: `τ ∝ 2√(R²-b²)` | ≤ 3,6×10⁻¹⁵ |
| T3 | b > R: τ = 0 exato, sem NaN | exato |
| T4 | halo p=2, forma fechada com e sem b | ≤ 4,9×10⁻¹⁶ |
| T5 | Terra: X ≈ 7,0×10⁹ g/cm², τ ≈ 4,2 | coerente |
| T6 | interpolação log-log exata em lei de potência; lança fora da faixa | 1,0×10⁻¹⁴ |

44 verificações, 0 falhas.

Para T4 a forma fechada com `b > 0` também foi derivada, e é o teste mais forte
do ponto de retorno:

```
τ = 2 N_A σ ρ₀ r₀² · (1/b)·[ arccos(b/r_out) − arccos(b/r_lo) ],  r_lo = max(b, r_in)
```

que no limite `b → 0` vira `2 N_A σ ρ₀ r₀² (1/r_in − 1/r_out)`.

## Ainda não implementado

Fase 2 (Schwarzschild, captura, redshift atuando) e Fase 3 (disco espesso,
varredura em paralelo, regeneração NC). Nada da Fase 2 está aqui: `metric` só
aceita `minkowski`, e `Metric::r_turning` tem uma implementação genérica por
bisecção que **não trata captura** — está documentado no ponto.
