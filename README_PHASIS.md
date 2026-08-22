# PHASIS — transporte DIS ao longo de geodésicas nulas

Profundidade óptica de neutrinos de ultra-alta energia ao longo de geodésicas
nulas em espaço-tempo estático e esfericamente simétrico, com perfil de matéria
arbitrário.

**Fases 1–2 (estas):** espaço plano e Schwarzschild. Sem disco, sem
dependência em θ, sem regeneração.

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

Implementações: `Minkowski`, `Schwarzschild`; `UniformBall`, `PowerLawHalo`;
`PowerLawCrossSection`, `TableCrossSection`, `ScaledCrossSection`,
`SumCrossSection`.

## Fase 2 — relatividade geral

**Ponto de retorno em forma fechada.** `b = r_t/√(1−r_s/r_t)` é a cúbica
`r_t³ − b²r_t + b²r_s = 0`, cujo discriminante é `b⁴(4b² − 27r_s²)`. A raiz
física é

```
r_t = (2b/√3)·cos[ (1/3)·acos( −3√3·r_s/(2b) ) ]
```

Exata, sem iteração. A bissecção genérica continua em `Metric::r_turning` como
fallback, e um teste confronta as duas (concordância ≤ 2×10⁻¹⁶).

**Captura.** O critério *é* o discriminante da cúbica: `b ≤ (3√3/2)r_s`. Não é
"o root-finder não convergiu" — `Metric::is_captured` tem override exato em
Schwarzschild. Retorna `Result{captured, τ=∞, P=0}`, nunca NaN.

**A conexão que fecha o argumento.** Avaliando `T` no próprio ponto de retorno:

```
T(r_t) = 2r_t − 3r_s
```

que zera **exatamente** na esfera de fótons `r_ph = 3r_s/2`. Não é coincidência:
ali `r_t` vira raiz **dupla** de `w`, e é por isso que `b_crit` existe. Duas
consequências, ambas exploradas no código:

1. O critério de captura e a estrutura da singularidade são o mesmo objeto.
2. Perto de `b_crit`, `w ~ (r−r_t)²` e o integrando vai como `1/s`. A
   substituição `s²` matou a raiz quadrada, **não** o logaritmo. `∫ds/s` diverge
   — e a divergência é **física**: o raio enrola na esfera de fótons e o
   comprimento próprio realmente cresce sem limite quando `b → b_crit⁺`.

**Tratamento do regime quase-crítico.** Com `ε = T(r_t)/r_t` medindo a distância
à criticalidade e `s_★ = √(T(r_t)/T'(r_t))` a largura do pico, a integral é
partida em `s_split = 8s_★`: abaixo dela, direto em `s`; acima, substituição
exponencial `s = s_split·e^u`, onde `ds/s = du` e a cauda `1/s` vira constante.
`Result.near_critical` sinaliza `ε ≤ 10⁻⁶`, e `winding_turns` reporta as voltas
além de uma linha reta. Nunca NaN, nunca truncamento silencioso.

## Corrente carregada ou total?

As tabelas de `dipole/` são **CC puras**. Sem regeneração, as duas convenções
são defensáveis por motivos **opostos**: `σ_CC` porque só a corrente carregada
remove o neutrino de vez; `σ_tot` porque, sem regeneração, a corrente neutra
também tira o neutrino do bin de energia.

Por isso `xsec_current` é explícito, vai para o cabeçalho do CSV, e **não tem
default silencioso**: `xsec_current=total` exige `table_path_nc` ou
`nc_to_cc_ratio` declarado.

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

## Testes de aceitação — Fase 2

| | o que verifica | resultado |
|---|---|---|
| T7 | desvio em relação a Minkowski escala como `(r_s/r_t)^1` | expoente 1,00006 |
| T8 | `Δφ → 2r_s/b`, erro caindo como `r_s/b` | erro ≤ 1,5×10⁻³ |
| T8c | coeficiente de 2ª ordem `= 15π/16` | 2,9455 vs 2,9452 |
| T9a | forma fechada vs bissecção genérica | ≤ 2,1×10⁻¹⁶ |
| T9b | `b_crit/r_s = 2,5980762114`; `r_t → 1,5r_s` como `√δ` | expoente 0,50007 |
| T9c | enrolamento diverge como `−ln δ`, sem NaN | inclinação `1/2π` |
| T10 | σ dentro vs fora da integral: diferença cresce como `(r_s/r_t)^1` | expoente 1,00005 |

38 verificações, 0 falhas. **T11** (não-regressão) é a suite da Fase 1 rodada
sem alteração: os mesmos 44 checks, nos mesmos dígitos.

Dois resultados saíram mais fortes do que o pedido e valem como validação
independente da geometria:

- **T8c.** O erro de T8 não é ruído: `Δφ = 2(r_s/b) + 2,9455(r_s/b)²`, e
  `15π/16 = 2,94524`. O integrador acerta o segundo coeficiente da expansão,
  não só o termo dominante.
- **T9c.** A inclinação medida das voltas contra `−ln δ` é `0,1591568`, e
  `1/2π = 0,1591549`. Isso significa `Δφ = −ln δ + const`, ou seja coeficiente
  de deflexão forte `ā = 1` — o valor analítico conhecido para Schwarzschild.
- **T9b.** `r_t/r_ph − 1 = 0,8165·√δ`, e `√(2/3) = 0,81650`: o expoente ½ e o
  coeficiente da raiz dupla, ambos.

Para T4 a forma fechada com `b > 0` também foi derivada, e é o teste mais forte
do ponto de retorno:

```
τ = 2 N_A σ ρ₀ r₀² · (1/b)·[ arccos(b/r_out) − arccos(b/r_lo) ],  r_lo = max(b, r_in)
```

que no limite `b → 0` vira `2 N_A σ ρ₀ r₀² (1/r_in − 1/r_out)`.

## Um achado no integrador, que vale registrar

`max_depth` limita a **profundidade** da recursão, não o **número de
avaliações** — e a recursão é binária, então 50 níveis são 2⁵⁰ folhas. Quando a
tolerância pedida cai abaixo do epsilon de máquina (`|S2−S1|` não consegue mais
encolher), o critério nunca é satisfeito e o processo trava. Apareceu na
primeira execução de T9c com `rel_tol = 10⁻¹⁴`.

`IntegratorOpts::max_evals` é o teto global que transforma isso em resultado
degradado e **sinalizado** (`Result::tolerance_met = false`) em vez de um
processo pendurado. Importa mais ainda quando forem 10⁵ geodésicas em paralelo.

## Ainda não implementado

Fase 3: disco espesso com `ρ(r,θ)`, varredura paralela, regeneração NC. O
caminho não esférico existe e lança `std::logic_error` explícito — a estrutura
dos dois ramos já está pronta, falta o ângulo `ψ` ao longo do raio.
