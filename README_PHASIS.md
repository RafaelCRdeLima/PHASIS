# PHASIS — transporte DIS ao longo de geodésicas nulas

Profundidade óptica de neutrinos de ultra-alta energia ao longo de geodésicas
nulas em espaço-tempo estático e esfericamente simétrico, com perfil de matéria
arbitrário.

**Fases 1–3 (estas):** espaço plano e Schwarzschild, disco com dependência em
θ, varredura paralela. Sem regeneração.

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

## Fase 3 — dependência em θ

**A geometria.** A métrica continua esfericamente simétrica, então o movimento
continua **planar** — basta girar o plano orbital. O plano tem normal a ângulo
`i` do eixo do disco, e `ψ` é a posição azimutal dentro dele, medida da linha de
nodos:

```
cos(θ) = sin(i)·sin(ψ)
```

`i = 0` dá `θ = π/2` sempre (raio no plano do disco); `i = π/2` cruza o disco em
`ψ = 0` e `π`. O raio é `(E_inf, b, i, ψ_t)`, com `ψ_t` no ponto de retorno.

**O acoplamento, e por que a arquitetura muda.** `ψ(r)` é necessário **dentro**
do integrando de τ, então as duas quadraturas deixam de ser independentes.
Sub-quadratura aninhada custaria O(N²) por raio — inviável para 10⁵. Com a mesma
fatoração das Fases 1–2, os dois integrandos ficam regulares em `s`:

```
dl/ds   = 2r     ·√( h·f(r_t)/T )
dψ/ds   = (2b/r) ·√( h·f·f(r_t)/T )
```

e o estado `(ψ, l, τ_in, τ_out, X_in, X_out)` avança junto num Dormand–Prince
5(4). Os dois ramos compartilham `ψ_off` e `l`; só o **sinal** com que `ψ_off`
entra em θ difere — `ψ_t − ψ_off` na entrada, `ψ_t + ψ_off` na saída. Nunca se
dobra um ramo.

Identidade de checagem cruzada, exata e verificada em T13:

```
dψ/dl = b·√f / r²
```

**Perfis:** `FlaredThinDisk` (`ρ = ρ₀(R/R₀)^{−p}·exp(−z²/2H(R)²)`,
`H = H₀(R/R₀)^q`, defaults p=15/8, q=9/8) e `QuasiSphericalADAF`
(`ρ ∝ r^{−3/2}`, ainda esférico, ponte entre os dois regimes).

**Paralelismo.** `sweep()` roda `#pragma omp parallel for` sobre um `vector<Ray>`
plano escrevendo em posições pré-alocadas. Sem RNG, sem estado mutável
compartilhado, sem I/O no laço — daí o determinismo bit a bit (T18). Duas coisas
obrigatórias: o corpo do laço fica em `try/catch` (exceção atravessando a
fronteira de OpenMP é UB), e a agregação é serial, depois do laço.

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

## Testes de aceitação — Fase 3

| | o que verifica | resultado |
|---|---|---|
| T12 | rota de EDO ≡ quadratura das Fases 1–2, perfil esférico | ≤ 1,7×10⁻¹³ |
| T13 | identidade `dψ/dl = b√f/r²` em 150 pontos | 3,6×10⁻¹⁶ |
| T14 | simetria `z→−z`: `τ(i,ψ_t) = τ(i,ψ_t+π)` | ≤ 6,8×10⁻¹⁵ |
| T15 | ramos genuinamente diferentes (8 configurações) | 50% a 200% |
| T16 | `i=0` ⇒ disco ≡ `PowerLawHalo`; limite `H→∞` | ≤ 4,6×10⁻¹³ |
| T17 | disco em campo plano, forma fechada | ≤ 6,7×10⁻¹³ |
| T18 | determinismo bit a bit com 1, 2, 4, 8 threads | exato |
| T18b | exceção não escapa da região paralela | 1000/1000 |

28 verificações, 0 falhas. **T19** (não-regressão): as suites das Fases 1 e 2
rodadas sem alteração, mesmos dígitos.

Uma nota sobre T16: o enunciado "com `H₀/R₀ = 10⁶` o disco converge para
`PowerLawHalo` com o mesmo p" só vale **no plano do disco**. Como `ρ` depende de
`R = r·sin θ`, o limite `H→∞` é `ρ₀(r·sinθ/R₀)^{−p}`, que coincide com
`PowerLawHalo` apenas em `θ = π/2`. O teste verifica as duas coisas
separadamente: a identidade exata em `i=0` (para qualquer H, o que é um teste
forte do mapeamento `(R,z)`) e a convergência para o limite correto em `i=0,9`.

## Dois achados no integrador, que valem registrar

`max_depth` limita a **profundidade** da recursão, não o **número de
avaliações** — e a recursão é binária, então 50 níveis são 2⁵⁰ folhas. Quando a
tolerância pedida cai abaixo do epsilon de máquina (`|S2−S1|` não consegue mais
encolher), o critério nunca é satisfeito e o processo trava. Apareceu na
primeira execução de T9c com `rel_tol = 10⁻¹⁴`.

`IntegratorOpts::max_evals` é o teto global que transforma isso em resultado
degradado e **sinalizado** (`Result::tolerance_met = false`) em vez de um
processo pendurado. Importa mais ainda quando forem 10⁵ geodésicas em paralelo.

**2. Controle de erro relativo puro colapsa o passo.** Na primeira execução de
T18, 180 de 10 000 raios saíam degradados. O sumário agregado é o que tornou
isso visível — com 10⁴ linhas ninguém percebe olhando o CSV.

Duas causas, ambas diagnosticadas medindo em vez de adivinhando:

- `abs_tol ≈ 0` no controlador da EDO. O ramo que se afasta do disco tem
  `τ ~ 10⁻⁸⁰`; com `sc = abs_tol + rel_tol·|y|` isso dá `sc ~ 10⁻⁸⁹`, e qualquer
  ruído nessa componente produz erro gigante. O passo encolhe até `h_min` por
  causa de uma quantidade fisicamente irrelevante. Um piso de `10⁻¹⁰` —
  desprezível para **todas** as componentes deste problema — eliminou 133 dos
  180.
- Os 47 restantes falhavam no trecho **vazio** `[0, s_lo]`, não onde há matéria.
  Os perfis têm corte duro em `r_in`, e o último estágio de Runge–Kutta desse
  trecho cai exatamente em `r = r_in`, onde `ρ` salta de 0 para finito. Degrau
  na fronteira: encolher o passo nunca ajuda. Como o trecho não tem matéria por
  construção, τ e X simplesmente não são integrados ali.

Resultado: 180 → 47 → **0**. Distinguir `hit_min_step` de `hit_max_steps` no
`OdeStats` foi o que permitiu separar os dois casos — as ações corretivas são
opostas.

## Um limite que vale saber

Todo raio não capturado tem `r_t ≥ r_ph = 3r_s/2`, e `f(3r_s/2) = 1/3`. Logo

```
E_loc^max / E_inf = 1/√f(r_t) ≤ √3
```

O blueshift em Schwarzschild é **limitado por √3, sempre**. Com `σ ∝ E^0.36` isso
é no máximo `3^0.18 ≈ 1,22` — 22% em σ, e só no ponto de retorno de um raio
quase-crítico. **A assinatura de RG no observável vem da geometria, não do
redshift:** o efeito grande é o comprimento de caminho dos raios quase-críticos,
que diverge logaritmicamente. (Em Kerr o limite muda; não é universal.)

Corolário prático: com `E_inf` dentro da faixa da tabela, `E_loc` nunca sai por
mais que √3. Meia década de folga no topo da tabela basta.

## Fase 4 — regeneração NC (parcial: etapas 1–3)

**O ponto que torna isso tratável.** Uma interação NC em raio `r` leva
`E'_loc → (1−y)E'_loc`. Como inicial e final estão no **mesmo** `r`, o mesmo
`√f(r)` aparece nos dois:

```
E_loc/E'_loc = E_inf/E'_inf = 1 − y
```

`y` é **invariante sob redshift**. O mapeamento entre bins de `E_inf` não depende
da posição; só o **coeficiente** varia ao longo do raio, porque σ é avaliada em
`E_loc`. Daí `CascadeKernel` montar os intervalos de `y` **uma vez** e guardar
`G_ij`, com `K_ij(E) = σ_nc(E)·G_ij` — uma multiplicação por passo, nenhuma
quadratura.

**Identidade de consistência**, exigida na construção (o construtor lança se
falhar): `Σ_i G_ij + G_leak_j = 1`. O vazamento é acumulado numa componente
própria do estado — nunca descartado, nunca empilhado no bin de baixo.

| | resultado |
|---|---|
| identidade discreta, `ConstantY` e `PowerLawY`, 5 a 40 bins/década | ≤ 8,2×10⁻¹⁵ |
| **T20** `dσ/dy = 0` ⇒ reproduz `exp(−τ)` das Fases 2–3 | 2,3×10⁻¹⁶ |
| **T21** conservação de número (`σ_CC = 0`) | 5,1×10⁻¹⁶ |
| **T22** solução analítica, perda **e** ganho juntos | ≤ 9,7×10⁻¹² |

**Sobre a forma de T22.** Com bins log-uniformes e fluxo em lei de potência, os
intervalos de `y` dependem só de `(i−j)`: a matriz é uma **convolução**, e a lei
de potência é **autovetor exato** do operador discreto, com

```
σ_eff,disc = σ_CC + σ_NC(1 − Z_disc),   Z_disc(γ) = Σ_d G_d ρ^{−d(γ−1)}
```

A solução da EDO é então exatamente `exp(−n σ_eff,disc L)` — é isso que fecha em
10⁻¹². Cobrar 10⁻¹⁰ contra o `σ_eff` **analítico** num número finito de bins
seria cobrar do esquema binado algo que ele não pode dar; a convergência
`Z_disc → Z = 1/γ` é medida à parte, e dá **ordem 2**:

```
  gamma  bins/dec   Z_disc          Z_analitico     desvio
  1.5    10         0.6681345510    0.6666666667    2.202e-03
  1.5    40         0.6667586944    0.6666666667    1.380e-04     ordem 1.99
  1.5    160        0.6666724195    0.6666666667    8.629e-06     ordem 2.00
```

### Piso do controlador: por componente, não global

Na cascata os bins de alta energia caem muitas ordens enquanto os de baixa
crescem, **no mesmo vetor de estado**. Um piso global escalado pelo maior
componente deixa os pequenos sem controle; escalado pelo menor, colapsa o passo.
A escolha foi um piso **por bin, proporcional ao fluxo inicial daquele bin**
(`10⁻¹⁴·φ_i(0)`): diz "não me importo com precisão relativa depois que este bin
caiu muito abaixo do que ele próprio começou" — que é fisicamente o certo, já
que um bin nessa situação não contribui mais.

### O invariante de borda, que custou dois bugs

`ρ` **tem de ser inclusiva nas duas bordas do suporte**: `ρ(r_support_min)` e
`ρ(r_support_max)` devem valer o limite pelo interior, nunca zero. Os
integradores usam exatamente `[r_lo, r_hi]`, então as bordas **são** pontos de
avaliação; se `ρ` for zero na borda e não-zero um fio para dentro, o estágio de
Runge–Kutta que cai ali vê um degrau, o estimador não converge, e encolher o
passo nunca ajuda porque a descontinuidade está no extremo.

Apareceu duas vezes, em sentidos opostos: na Fase 3 um trecho **vazio** tocava
`r_in` (onde ρ é corretamente não-nula) — resolvido não integrando matéria onde
não há; na Fase 4 o ramo de entrada **começava** em `r_out`, onde `UniformBall`
devolvia zero por usar `r < R` em vez de `r <= R`. Está documentado como
invariante em `DensityProfile`.

## Ainda não implementado

Fase 4, etapas 4–5: **T23** (exponencial de matriz como rota independente),
**T24** (convergência de grade com expoente reportado), **T25** (magnitude do
redshift), round-trip do leitor de tabela, e a saída CSV do espectro.

O `TableDifferentialCrossSection` ainda não existe: as tabelas de `dσ_NC/dy` do
`dipole` **também não**. Ver a nota de bloqueio abaixo.

### Bloqueio conhecido: dσ_NC/dy não existe no `dipole`

`makeNCParameters` é chamada num único lugar do `dipole` (`main.cpp`, modo
`wavefunctions`, para desenhar `|ψ|²`). **Nada calcula σ_NC**, e só existe
`d2sigma_dxdy_CC`. Produzir a tabela exige: corrigir D10 (o α do NC está 2×
baixo), implementar as funções de estrutura NC, o duplo-diferencial, e a saída em
grade `(E, y)` com `y` logarítmico.

Isso **não bloqueia** as etapas 4–5: T23–T25 usam kernels analíticos, exatamente
como T20–T22.
