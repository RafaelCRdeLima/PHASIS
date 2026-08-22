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

### Etapas 4–5

| | resultado |
|---|---|
| borda: invariante genérico varrendo o **registro** de perfis | 4/4 |
| **T23** exp(nLM) vs EDO, condição inicial de **bin único** | ≤ 1,0×10⁻¹⁵ |
| **T23b** cada `K_ij` lida pela expansão de L curto | 3,4×10⁻¹⁰ |
| **T24** ordem observada por janela, `ConstantY` e `PowerLawY` | ver abaixo |
| **T25** redshift: `M(E_loc)` vs `M` congelada, mesma geometria | 4,1% em τ |
| round-trip do leitor + metadados obrigatórios | 1,8×10⁻¹⁴ |

**T23 usa bin único, não lei de potência**, e a razão é estrutural: a lei de
potência é autovetor **exato** do operador discreto, então T22 testa **uma**
direção num espaço de N dimensões. Um erro que preserve `σ_eff` para leis de
potência mas erre a distribuição entre bins passaria limpo. T23b vai além e lê
**cada entrada** `K_ij` — a identidade de consistência só testa `Σ_i K_ij`.

**T24 — a transição prevista aparece onde deve.**

```
    ConstantY  (suave)
    n/dec    ln(rho)     erro rel        ordem local
    5        0.4605      2.5936e-02         -
    40       0.0576      4.1406e-04       1.999
    640      0.0036      1.6080e-06       2.007        ordem media 1.997

    PowerLawY  beta = 1,  y_min = 0.05
    n/dec    ln(rho)     erro rel        ordem local  regime
    5        0.4605      5.6231e-02         -         sub-resolvido
    20       0.1151      4.4068e-03       1.856       sub-resolvido
    40       0.0576      9.8379e-04       2.163       sub-resolvido
    80       0.0288      1.3147e-04       2.904       corte resolvido
    640      0.0036      1.5982e-06       3.474       corte resolvido

      ordem sub-resolvida (n/dec  5 ->  20) : 1.837
      ordem resolvida     (n/dec 160-> 640) : 2.915
      transicao prevista em ln(rho) = y_min = 0.050  =>  n/dec = 46
```

A grade em `y` induzida é `y_k = 1 − ρ^{−k}`, com espaçamento `≈ ln ρ` perto de
zero. A transição cai entre n/déc 40 (`ln ρ = 0,0576 > y_min`) e 80
(`ln ρ = 0,0288 < y_min`) — exatamente onde `ln ρ` cruza `y_min`. Ajustar uma
única reta global daria um expoente intermediário sem significado.

*Ressalva honesta:* no regime resolvido a ordem **local** oscila (0,53; 2,36;
3,47). O ajuste em janela larga dá 2,92, acima de 2 — sinal de que ali duas
contribuições de erro com sinais opostos se cancelam parcialmente. O que o teste
cobra é a **melhora** ao cruzar `y_min`, não um valor específico.

**T25 mede τ, não P** — e isso qualifica uma afirmação anterior deste README:

```
  b/r_s   E_loc/E_inf   tau      dif em tau   dif em P
  3.0     1.34730       9.654    0.0406       0.3239
  10.0    1.05747       3.749    0.0131       0.0479
  100.0   1.00506       0.340    0.0012       0.0004
```

O limite `√3` controla **σ**: a diferença em τ é 4,1% no pior caso, dentro do
teto absoluto `3^0.18 − 1 = 21,9%`. Mas **P não herda esse teto**:

```
δP/P ≈ δτ = τ · (δσ/σ) ≤ 0,219 τ
```

Na janela observável `τ ∈ [0,3; 5]` isso vai de ~7% a ~67%. Ou seja: **o redshift
não é desprezível exatamente onde o observável vive** — a amplificação é maior
justamente onde τ é grande o bastante para importar e pequeno o bastante para P
não ser zero. Uma versão anterior deste README dizia que "a assinatura de RG vem
da geometria, não do redshift"; isso está **errado** na janela que interessa.

### Uma predição testável, para quando a tabela do dipole existir

O teto do redshift sobre σ é `3^{α/2}`, com `α = d ln σ/d ln E`. Mas α **não é
constante**: redshift e saturação se acoplam com sinais opostos.

```
blueshift  ⇒  E_loc ↑  ⇒  x_loc ↓  ⇒  mais fundo na saturação  ⇒  supressão ↑
```

Com `α = 0,36` (colinear) o teto é 1,219; com `α ≈ 0,30` (saturado) cai para
1,176. **A própria saturação reduz o realce por redshift.** Um cálculo colinear
com α fixo superestima o efeito de RG; um cálculo de dipolo não, porque α emerge
da dinâmica.

Isso é mensurável neste código assim que a tabela existir: extrair `α_eff` da
tabela e verificar se o teto observado em T25 cai **abaixo** de 1,219. Se cair,
é um acoplamento GR×saturação demonstrado.

O código já tem o que precisa — `freeze_redshift` isola o efeito com a geometria
fixa, e T25 mede a diferença em τ. Falta só a tabela.

**Metadados obrigatórios no formato de tabela.** As nove chaves
(`convention_y`, `target`, `projectile`, `current`, `units_sigma`, `units_E`,
`M_Z_GeV`, `dipole_model`, `generated_by`) são exigidas: remover qualquer uma faz
o leitor lançar (9/9 verificadas), e `convention_y`/`current` são validadas
contra o que o chamador declara esperar. Não é burocracia — nada neste projeto
força a explicitar a convenção, porque quem escreve e quem lê é o mesmo autor, e
a convenção é onde mora o fator 2.

### O invariante de borda virou teste, não comentário

Duas ocorrências em fases diferentes é assinatura estrutural. O registro de
perfis é **auto-registrável** (`PHASIS_REGISTER_PROFILE`), e o teste varre o
registro comparando `ρ` na borda com `ρ` um passo para dentro. Perfil novo entra
no teste sozinho — inclusive os da Fase 5.

## Fase 5 — infraestrutura de varredura (parcial)

Com ~25 ordens de grandeza de variação em τ, grade uniforme gasta quase tudo
onde a resposta é 0 ou 1. A física vive na casca fina onde τ ~ 1 — e é só ali
que existe sensibilidade ao modelo de σ. `locate_shell` trabalha em `log τ`
sobre um parâmetro escalar qualquer (i, b, o que for).

**A armadilha que estrutura o código:** bissecção assume monotonicidade, e
`τ(parâmetro)` **não é monotônico**. Um perfil oco já basta — a corda por uma
casca `[r₁,r₂]` **cresce** com `b` até `b = r₁` (de `2(r₂−r₁)` a
`2√(r₂²−r₁²)`) e só depois cai. Bissecção cega devolveria uma raiz em silêncio,
e esse é o único modo de falha que não produz erro visível.

Protocolo, nesta ordem: varredura grosseira (≥32 pontos) contando mudanças de
sinal → exatamente uma, bissecta no bracket → mais de uma, registra **todas** e
marca `multi_root` → nenhuma, `shell_found = false` com o valor extremo, sem
extrapolar e sem NaN.

| | resultado |
|---|---|
| **T27** casca analítica, `b* = √(R² − (2N_Aσρ)⁻²)` | 3,2×10⁻¹³ |
| **T28** perfil oco: duas raízes, ambas achadas, `multi_root` marcado | ≤ 2,1×10⁻¹³ |
| **T29** sem bracket (transparente / opaco), retorno limpo | sem NaN |

Em T28 o localizador acha `b = 8,182×10⁷` e `1,561×10⁸`, com τ = 1 nas duas a
10⁻¹³, e sinaliza a não-monotonicidade.

## Ainda não implementado

- **O ensemble não está decidido** — ver abaixo. Nada da infraestrutura acima
  depende dele.
- T30 (determinismo da varredura completa), cabeçalho de reprodutibilidade
  (commit, hashes das tabelas, `OMP_NUM_THREADS`), checkpointing.
- Saída CSV do espectro transmitido.

## Topologia do raio emitido (T31–T35)

**Decidido: emissão.** Emissor **estático** em `r_emit`; emissor orbitando é a
fase seguinte.

O mapeamento ângulo → `b`, com `ψ` medido da direção radial para fora:

```
b = r_emit · sin ψ / √f(r_emit)
```

Emissão isotrópica no referencial local é **cos ψ uniforme** em [−1,1], não ψ
uniforme.

### São quatro casos, não três

`g(r) = r/√f(r)` tem `g'(r) = ½ r^{1/2}(r−r_s)^{−3/2}(2r − 3r_s)`, que zera na
esfera de fótons. **Abaixo dela `g` decresce**, então um raio emitido para fora
encontra `g(r) = b` *acima* de `r_emit`, vira, e cai. A condição de escape
inverte:

| `r_emit` | direção | escapa se |
|---|---|---|
| `> 1,5 r_s` | fora | sempre |
| `> 1,5 r_s` | dentro | `b > b_crit` |
| `< 1,5 r_s` | fora | **`b < b_crit`** ← invertida |
| `< 1,5 r_s` | dentro | nunca |

Em `r_emit = 1,5 r_s` os dois ramos concordam (ali `b ≤ b_crit` sempre), então
não há descontinuidade — e é exatamente onde `data/exemplo_buraco_negro.cfg`
coloca `r_in`.

### O cone de escape, com dois valores exatos

```
sin ψ_c = b_crit·√f(r)/r        f_cap = (1 − cos ψ_c)/2
```

| | resultado |
|---|---|
| `ψ_c(1,5 r_s) = π/2`, `f_cap = 1/2` | rel 0 |
| `ψ_c(3 r_s) = π/4`, `f_cap = (2−√2)/4` | rel 0 |
| fronteira de `classify()` vs `π − ψ_c`, 6 raios | 2,2×10⁻¹⁶ |
| Monte Carlo, 10⁶ amostras, 6 raios | pior 0,01 σ |
| **T33** emitido em `r_out` ≡ vindo do infinito | **0** |
| **T34** radial `b=0` vs `b` pequeno | ≤ 4,9×10⁻¹⁶ |

Na esfera de fótons, **exatamente metade** do céu local é captura. Na ISCO,
**exatamente 45°** de cone e 14,64%.

### Um caso que a especificação não previa

Para `OutboundOnly` com `b < b_crit` o ponto de retorno **não existe** — não é um
ponto de retorno "virtual" que o raio não alcança, é ausência. Mas aí `w = 1 −
f b²/r²` nunca se anula na faixa percorrida, então **não há singularidade e não
há o que fatorar**: integra-se direto. São três modos:

| modo | quando | `dl/ds` |
|---|---|---|
| `Normal` | há retorno e o raio o alcança | `2r√(h·f(r_t)/T)` |
| `NoTurn` | emitido para fora, sem retorno | `2s√(h/w)`, âncora em `r_emit` |
| `Radial` | `b = 0` | `2s√h` — `r_t` e `f(r_t)` nem são invocados |

O caso radial precisa de caminho próprio porque em Schwarzschild `r_t = 0` e
`f(0)` são ambos singulares. Também foi corrigido um guarda `b > 0` que fazia um
raio **radial vindo do infinito** escapar da checagem de captura — em
Schwarzschild ele cai no buraco negro.

### O emissor orbitando não é uma correção pequena

`β = √(M/r)/√(1−2M/r)`, e na ISCO (`r = 6M`) isso vale **exatamente 1/2**, com
`γ = 2/√3`. O fator Doppler `γ(1 ± β)` varre `1/√3` a `√3` — **fator 3 em
energia**, contra o teto `√3` do redshift gravitacional. A aberração do fluido é
o **maior** dos dois efeitos. Adiá-la é sequenciamento, não um argumento de que
seja pequena.

(Uma versão anterior deste README estimou `β ≈ 0,35`. A fórmula estava certa mas
foi avaliada em `r = 10M`, não na ISCO.)

### A decisão que estava pendente: emissão ou transmissão

**(a) Emissão** — fonte em `r_emit` dentro do disco, isotrópica no referencial
local. É a física: o AGN produz os neutrinos lá dentro. Mas quebra a topologia
atual: aparecem raios de **ramo único** (emitidos para fora), e um caso novo
(emitido para dentro com `b < b_crit`, que cai no buraco negro).

O mapeamento ângulo local → `b` para um emissor **estático** é limpo:

```
b = r_emit · sin ψ / √f(r_emit)
```

com `ψ` medido da direção radial. Note que `b_max = r_emit/√f(r_emit)` em
`ψ = π/2` é exatamente a relação do ponto de retorno — um raio emitido
tangencialmente tem seu ponto de retorno em `r_emit`. Para um emissor
**orbitando** há aberração adicional entre o referencial do fluido e o estático,
que é uma decisão à parte.

**(b) Transmissão** — de ∞ a ∞, sombra do disco. Continua funcionando sem
mudança: `r_emit_cm ≤ 0` é o default e significa "vindo do infinito".

**Escolhido (a).** T33 garante que a máquina antiga não quebrou: um raio emitido
em `r_out` reproduz o resultado de ∞ a ∞ para o mesmo `b` **exatamente** (dif
0,00e+00), e emitido no meio do halo tem o ramo de saída idêntico e o de entrada
truncado.

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

## Acoplamento GR × saturação (T36–T39)

A grandeza é propriedade **só de `σ_tot`**. Não depende de NC, nem de `dσ/dy`,
nem do acoplamento NC do `dipole` (D10). Como `σ_CC` já estava validada, a
medida estava desbloqueada.

```
α_eff(E) = d ln σ / d ln E        teto(E) = 3^(α_eff(E)/2)
```

O teto vem de `E_loc/E_inf = 1/√f(r_t) ≤ √3` em Schwarzschild: todo raio não
capturado tem `r_t ≥ 1.5 r_s`, onde `f = 1/3`.

**A assinatura não é o valor do teto, é a derivada dele em E.** Um deslocamento
de normalização qualquer reajuste de PDF reproduz; uma inclinação que cai com a
energia, não.

### O estimador — e por que não é analítico

`CrossSection::log_slope_detail` usa **diferença centrada em ln E** com passo
adaptativo por redução pela metade (h inicial 0,1; piso 3e-3).

A rota analítica está **fechada por construção**: `TableCrossSection` interpola
linearmente em (ln E, ln σ), o que é C⁰. A derivada analítica dessa spline é uma
função escada, com salto em cada nó, e no próprio nó nem está definida.

O critério de parada é "achei o segmento", não "reduzi o truncamento": enquanto
a janela `[E e^-h, E e^+h]` contém um nó, `D(h)` é média ponderada de dois
segmentos e muda ao encolher; assim que a janela cabe dentro de um segmento,
`D(h)` para de mudar porque ali `ln σ` é exatamente linear em `ln E`.

O denominador é `ln(E₊) − ln(E₋)`, calculado nos mesmos valores que `σ` recebeu,
não `2h`. `E·e^h` seguido de `ln` não devolve `h` exatamente, e usar `2h`
injetaria erro relativo `~1e-16/h` — em `h = 3e-3` isso é 3e-14, na mesma ordem
do que T36 mede.

`log_slope` **não é virtual**, de propósito: se `PowerLawCrossSection`
devolvesse `α` de um campo, T36 — o único teste que valida o estimador — estaria
validando um `return`.

### T37 — uma identidade, não um ajuste

`R(b,E) = (τ_full − τ_frozen)/τ_frozen = ⟨f^(−α/2)⟩ − 1`, média com peso
`n(r)σ dl`. Daí, com `g = f^(−α/2)` e `g_max = 3^(α/2)`:

```
teto − 1 − R  =  ∫w(g_max − g) / ∫w  =  K / τ_frozen
```

O numerador **converge** quando `b → b_crit⁺`, porque `(g_max − g)` se anula
exatamente onde o peso diverge; o denominador diverge logaritmicamente. Então
dois pontos bastam para extrapolar, resolvendo duas equações numa reta:

```
teto − 1 = (R₁τ₁ − R₂τ₂)/(τ₁ − τ₂)
```

Medido, com `K` constante a 6 algarismos ao longo de 6 décadas em `b/b_crit − 1`:

| α | K | 1+R extrapolado | 3^(α/2) | erro rel |
|---|---|---|---|---|
| 0,20 | 8,55704 | 1,116123174 | 1,116123174 | 3,0e-11 |
| 0,30 | 21,18490 | 1,179147646 | 1,179147646 | 4,6e-11 |
| 0,36 | 34,34048 | 1,218657950 | 1,218657950 | 5,5e-11 |

**Correção aos valores da especificação:** os tetos são **1,11612**, **1,17915**
e **1,21866**. A especificação dizia 1,128 e 1,176 para os dois primeiros;
`3^0,10 = 1,11612` e `3^0,15 = 1,17915`. O terceiro (1,219) estava certo.

O halo do teste tem `r_in = 1,2 r_s`, **abaixo** da esfera de fótons: sem isso os
raios quase-críticos enrolariam no vazio e o teto seria inatingível.

### T39 — o pré-requisito lógico

`R` é invariante sob `σ → kσ` porque um fator global multiplica numerador e
denominador. Medido: **exatamente 0** para `k ∈ {0,5, 2}` (potências de dois, a
reescala é exata em binário) e `≤ 2,6e-16` para `k = 10`; pior caso de 27
combinações (3 k × 3 b × 3 E) igual a **2,1e-13**, e 4,7e-15 na tabela real.

Isso vale porque `adaptive_simpson` usa `abs_tol = 0`: o critério é puramente
relativo, então escalar o integrando escala `S1`, a tolerância e `|S2−S1|` na
mesma proporção e o **padrão de subdivisão é idêntico**.

É isso que legitima comparar dipolo com colinear apesar da supressão global do
primeiro: `R` isola a **forma**.

### T38 destapou um bug no `dipole`

Rodado na tabela de produção de então, o teto pontual **não** era monótono: 92
segmentos subiam de 189. O ruído de `ln σ` na tabela era **8,2e-4**, e ele
escalava como **N⁻¹** com a densidade de nós da quadratura — não como N⁻⁴.
Duplicar a tabela de `F` de 121 para 241 nós não mudava **nada**: as duas
corridas saíram idênticas byte a byte.

`N⁻¹` com a tabela de `F` irrelevante só tem uma causa: **descontinuidade no
extremo da quadratura**. Em `sigma_nuN_core.cpp`,

```cpp
if (y <= 0.0 || y >= 1.0) return 0.0;   // ERRADO em y = 1
```

O laço em `x` integra a partir de `x = Q²/s`, que é **exatamente `y = 1`**. Mas
`y = 1` é o extremo cinemático (lépton de saída com energia nula), não um ponto
proibido: ali o integrando vale `prefator·[F₂/2 − F_L/2 + xF₃/2]`, finito e da
ordem do resto. Zerá-lo punha um degrau no nó de borda, e Simpson com o valor da
borda errado erra em O(h), não em O(h⁴).

É a **mesma classe** dos dois bugs de borda de densidade das Fases 3 e 4:
fronteira dura avaliada exatamente num extremo de integração. Terceira
ocorrência.

Depois do conserto, medido em 80 energias entre 1e7 e 1e14 GeV:

| | ruído de `ln σ` | subidas de α |
|---|---|---|
| antes, 16 nós/década | 8,2e-4 | 37 de 78 |
| antes, 48 nós/década | 2,6e-4 | 26 de 78 |
| **depois, 16 nós/década** | **1,2e-6** | **0 de 78** |
| **depois, 48 nós/década** | **1,2e-6** | **0 de 78** |

Fator **670** a 16 nós/década. E 16 contra 48 nós/década agora diferem em no
máximo **5,2e-7**: a quadratura está convergida, e o 1,2e-6 residual é do meu
polinômio de referência, não do integrador.

O conserto também remove um **viés**: `σ` sobe 0,23% em média (até 0,39%) a 16
nós/década, 0,076% a 48. O viés era o mesmo termo de borda faltando; o ruído era
a parte dele que oscilava conforme `Nlogx` saltava de 2 em 2 com a energia.

A tabela de produção foi regerada. Nas duas (GBW e bCGC), **0 de 189** segmentos
sobem em `E ≥ 1e7 GeV`.

| E (GeV) | α GBW | teto GBW | α bCGC | teto bCGC |
|---|---|---|---|---|
| 1e5 | 0,6346 | 1,4171 | 0,6341 | 1,4167 |
| 1e7 | 0,3565 | 1,2164 | 0,4112 | 1,2534 |
| 1e9 | 0,2796 | 1,1660 | 0,3008 | 1,1797 |
| 1e11 | 0,2561 | 1,1511 | 0,2434 | 1,1431 |
| 1e13 | 0,2439 | 1,1434 | 0,2085 | 1,1214 |

O teto cai de 1,22 para 1,14 (GBW) e de 1,25 para 1,12 (bCGC) entre 1e7 e 1e13.
Essa queda é a assinatura.

### Validação cruzada de metadados

Toda tabela de `σ_tot` carrega as **nove chaves obrigatórias**, no formato
`# chave = valor`, mais o hash do commit que a gerou (com sufixo `-sujo` se a
árvore estava modificada).

`phasis::assert_comparable(a, b)` exige as nove nas duas e **igualdade** de
`target`, `projectile`, `current`, `units_sigma` e `units_E`. Não exige
igualdade de `dipole_model` nem de `generated_by` — é exatamente ali que as duas
devem diferir.

Este é o único ponto onde uma comparação errada passaria despercebida: cada
tabela isolada é autoconsistente, e a razão entre uma de ν e uma de ν̄ continua
sendo um número perfeitamente calculável. Só a comparação é que deixa de
significar algo, e nada **dentro** de cada arquivo consegue detectar isso.

### A tabela colinear

Gerada fora do C++, com **yadism 0.13.11 + eko 0.15.5** (NLO, ZM-VFNS, TMC
desligada, CKM completa) sobre **NNPDF3.1 NLO**, tabulando `F₂`, `F_L` e `xF₃`
numa grade de 121 × 71 nós em `(x, Q²)` — `x` de 1e-15 a 1, `Q²` de 1 a 2e14.

Convenção conferida contra a fórmula de ordem dominante antes de gastar CPU: em
yadism o observável chamado **`F3` é xF₃**, não F₃. A ordem dominante do yadism
reproduz `2x[d+s+b−ū−c̄]` no nó, e `F_L = 0` exatamente.

A integração em `σ` usa a **mesma fórmula mestra** do `dipole`
(`tools/make_collinear_table.py` reproduz `sigma_nuN_core.cpp` linha por linha,
inclusive a regra de nós por década e o `Nlogx` fixo ao longo do laço de `Q²`).
Duas diferenças deliberadas, e só duas: as funções de estrutura são colineares, e
**não há fator `(1−x)⁷`** — essa é a prescrição de grande `x` do modelo de
dipolo, e os PDFs colineares já se anulam em `x → 1` sozinhos.

Confere com a literatura colinear: `σ_CC(1e6 GeV) ≈ 6e-34 cm²`.

**A extrapolação está registrada, não escondida.** Abaixo de `x = 1e-9` o LHAPDF
extrapola. A tabela carrega uma quarta coluna com a fração de `σ` que vem dali:

| E (GeV) | fração de σ com x < 1e-9 |
|---|---|
| ≤ 2,3e10 | 0,0% |
| 1,9e11 | 0,5% |
| 1,6e12 | 4,3% |
| 1,3e13 | 20,7% |

Acima de ~1e12 GeV, `α_eff` colinear é cada vez mais uma **propriedade da
extrapolação** (uma lei de potência continuada), não uma medida. É exatamente a
hipótese que o dipolo substitui — mas quem ler o arquivo daqui a seis meses
precisa saber onde ela começa a mandar.

### O resultado

`teto_colinear − teto_dipolo`, em pontos percentuais de `σ`:

| E (GeV) | GBW | bCGC |
|---|---|---|
| 1e5 | **−6,5** | −6,4 |
| 1e7 | +3,4 | −0,4 |
| 1e9 | +3,6 | +2,2 |
| 1e11 | +2,3 | +3,1 |
| 1e13 | +1,1 | +3,3 |

O sinal **inverte** perto de 1e6 GeV: abaixo dali o dipolo sobe mais rápido que o
colinear e o teto dele é maior. A separação no sentido previsto — teto do dipolo
**abaixo** do colinear — só se estabelece na região de saturação, e vale 2 a 4
pontos percentuais, não os 4,3 estimados a priori.

`R(b,E)` medido diretamente, em `E = 1,09e7 GeV`, com `r_in = 1,2 r_s`:

| b/b_crit − 1 | r_t/r_s | R dipolo | R colinear |
|---|---|---|---|
| 1e-12 | 1,50000 | 0,1914 | 0,2237 |
| 1e-6 | 1,50123 | 0,1750 | 0,2044 |
| 4e-3 | 1,58275 | 0,1337 | 0,1555 |
| 1,0 | 4,59627 | 0,0281 | 0,0323 |

Consistente com T37: `3^(α/2) − 1` vale 0,2154 (dipolo, `α = 0,354`) e 0,2469
(colinear, `α = 0,404`), e o `R` medido se aproxima de cada um por baixo, como a
lei `K/τ_frozen` exige.

**A resposta honesta sobre detectabilidade** ainda não foi escrita, e ela é uma
pergunta separada: 3 pontos percentuais em `σ` viram `Δτ = 0,15` em `τ = 5`, ou
16% em `P`. Isso é detectável no cálculo. Se é detectável num telescópio depende
do ensemble, que ainda não está decidido.

### Executáveis

```bash
make                      # inclui alpha_scan e redshift_scan
./build/alpha_scan    --dipole A.dat --colinear B.dat --out data/alpha_scan.csv
./build/redshift_scan --dipole A.dat --colinear B.dat --out data/redshift_scan.csv
```

Os dois chamam `assert_comparable` antes de qualquer conta, e ecoam as nove
chaves das **duas** tabelas no cabeçalho do CSV.

Para regerar a tabela colinear (precisa do env micromamba `dis`):

```bash
export PATH=/home/rafael/micromamba/envs/dis/bin:$PATH
export LD_LIBRARY_PATH=/home/rafael/micromamba/envs/dis/lib:$LD_LIBRARY_PATH
python tools/build_collinear_sf.py      # ~13 min, RSS constante em ~800 MB
python tools/make_collinear_table.py    # ~5 s
```

O `build_collinear_sf.py` fatia o cálculo por `Q²` e grava cada fatia assim que
sai. Não é elegância: `yadism.Runner.get_result()` guarda a grade de funções de
coeficiente de **cada** ponto pedido até o fim da chamada — 266 MB + 0,95 MB por
ponto. Pedir os 8591 pontos de uma vez custa 8,4 GB, e numa máquina de 15 GB com
2 GB de swap isso trava o sistema inteiro sem disparar o OOM-killer. Aconteceu.
Fatiado, o consumo é constante e uma interrupção custa uma fatia.

## Corrente neutra (T44–T50)

Fecha o **P6**, e resolve o **P3** no caminho.

### D10 — duas metades que pareciam se contradizer

O comentário antigo dizia, na mesma respiração, que `α_NC` estava um fator 2
baixo e que a checagem numérica dava o valor certo do SM. As duas coisas são
verdade, e são **fatores 2 diferentes**.

**O valor estava errado.** Escrevendo o vértice como `C γ^μ(g_V − g_A γ₅)`, o Z
dá `C = g_Z/2 = g/(2 cos θ_W)`, e a convenção deste código é `α_EW = C²/4π`. Com
`g² = 4√2 G_F M_W²` e `M_W = M_Z cos θ_W`:

```
α_NC = √2 G_F M_Z² / 4π
```

O antigo era `G_F M_Z²/(√2·4π)` — exatamente metade. Confere com o CC pela mesma
rota: `C_CC = g/(2√2)` dá `α_CC = G_F M_W²/(√2·4π)`, que é o que estava escrito e
está certo.

**E não muda σ.** `StructureTable::computeAt` divide `F_T` e `F_L` pelo `α_EW`
**do próprio canal** (`structure_table.cpp:140`) — é a convenção DIS: toda a
normalização eletrofraca de σ vem de `G_F²` e do propagador na fórmula mestra, e
o que sobra em `F` é a soma efetiva de cargas. Então `α_NC` cancela.

Por isso a checagem antiga estava certa mesmo com o valor errado: ela mede a
razão das **somas de carga**, onde `α` já saiu.

### Uma fórmula mestra, duas correntes

`d2sigma_dxdQ2_param` e `d2sigma_dxdy_param` chamam o mesmo núcleo. A única
constante que distingue CC de NC é `M_V` (`bosonMass2`) e de onde vem `xF₃`.
Escrever as duas como uma função não é economia de linhas: é a garantia de que
não possam divergir — a razão `σ_NC/σ_CC` só significa algo se o resto for
idêntico, e aqui é o mesmo código.

O refactor foi verificado: a tabela CC saiu **bit a bit idêntica** à anterior, em
300 energias.

O CC junta dois sabores num dipolo (`ud̄`, `cs̄`): 2 canais. O NC junta o mesmo
sabor consigo (`uū`, `dd̄`, `ss̄`, `cc̄`): 4 canais. Mesma população de quarks,
contagem diferente de dipolos — e é daí que vem a maior parte da razão, não do
acoplamento.

```
Σ_NC (g_V² + g_A²) sobre u,d,s,c = 1,3127     (= o C de KK eq. 15)
Σ_CC (g_V² + g_A²) sobre ud, cs  = 4
razão × (M_Z/M_W)² = (1,3127/4)(1,2870) = 0,4224
```

Medido nas tabelas de produção:

| E (GeV) | σ_CC (cm²) | σ_NC (cm²) | razão |
|---|---|---|---|
| 1e3 | 3,741e-36 | 1,140e-36 | 0,3047 |
| 1e5 | 2,047e-34 | 7,335e-35 | 0,3584 |
| 1e7 | 1,777e-33 | 7,025e-34 | 0,3953 |
| 1e9 | 7,400e-33 | 3,001e-33 | 0,4055 |
| 1e11 | 2,521e-32 | 1,032e-32 | 0,4093 |
| 1e13 | 7,961e-32 | 3,281e-32 | 0,4122 |

A razão **cresce** com a energia porque o `xF₃` do CC domina embaixo (é de
valência: +70% em 1e3 GeV, +0,16% em 1e14) e some em cima. O assintótico,
0,409, fica 3% abaixo dos 0,4224 da contagem de cargas — a diferença é efeito de
massa do charme, que a contagem ignora.

### Um bug meu, e o que o pegou

A primeira tabela NC deu `σ_NC/σ_CC = 2,49` em vez de 0,42. Os canais eram
`uu, dd, ss, cc` mas os **acoplamentos continuaram CC** (`g_V=−1, g_A=+1`),
porque o argumento `current` simplesmente não foi passado ao construtor da
tabela — e ele tinha valor padrão `CC`.

O resultado era plausível: tabela lisa, monotônica, sem exceção, com `α_eff`
bem-comportado. Só o valor era errado. O que denunciou foi o **teste de física**
— σ_NC/σ_CC tem de dar 0,42 — e não nada interno ao código.

O conserto estrutural foi tirar o valor padrão. O compilador imediatamente
apontou o segundo call site (`test_table.cpp:55`), que eu não teria lembrado de
olhar. Um default silencioso num parâmetro que muda a **física** é um convite.

### A tabela dσ/dy, e o corte cinemático

`dσ/dy` só é não nula para `y ≥ Q2min/(s·x_max)`, e esse limite **depende da
energia** (cai como 1/E). O formato de tabela exige uma grade em `y` comum a
todas as energias.

A escolha: `y_min` da grade é o corte na energia **mais alta**, e nas energias
mais baixas as linhas abaixo do próprio corte saem zeradas. São zeros de verdade
— consequência de `Q² ≥ Q2min` — não lacunas. A alternativa esconderia, nas
energias altas, a faixa de `y` pequeno que a cascata mais usa: `y` pequeno é
perda de energia pequena, que é o regime onde o kernel de regeneração pesa.

O cabeçalho reporta, por energia, que fração de σ a integral em `y` recupera. Com
120 nós log-espaçados: 0,9988 em 1e3 GeV, **1,0104** em 1e14. O 1% é o erro do
trapézio, e está medido, não estimado.

### A "Fase 4b": kernel com forma dependente de E

`CascadeKernel` **lançava** para qualquer seção de choque com forma em `y`
dependente de E. A tabela real depende — e muito:

```
g(y) = (dσ/dy)/σ_NC        E=1e5      E=1e8      E=1e11
   y = 1e-3                1,326     43,342     53,621
   y = 1e-2                1,897     10,414      9,506
   y = 0,5                 0,843      0,373      0,308
```

Fator **40** de variação entre 1e5 e 1e11 GeV. Congelar a forma seria erro
grosseiro.

A rota construída: tabelar `G_ij` em nós **log-E** e interpolar **linearmente**
em `ln E`. A linearidade não é conveniência:

> A identidade de consistência `Σ_i G_ij + G_leak_j = 1` é **linear** nos `G`.
> Interpolação linear de um conjunto que soma 1 em cada nó ainda soma 1 entre os
> nós.

Ou seja, a cascata conserva número **exatamente em qualquer E**, não só nos nós.
Com interpolação cúbica isso deixaria de valer, e o erro apareceria como fluxo
criado ou destruído do nada ao longo do raio — o tipo de erro que nenhuma
checagem a jusante pegaria.

Medido: resíduo **2,22e-16** nos nós *e entre* eles. Seis nós em `ln E` já
resolvem a forma (diferença de 3,9e-3 contra doze nós).

Os intervalos em `y` que ligam o bin `j` ao bin `i` **não** dependem de `E_loc`:
`y` é invariante sob redshift e as bordas dos bins estão em `E_inf`. Só o
argumento de `dσ/dy` varia. É por isso que dá para tabelar em E de uma vez e
nunca mais reintegrar.

**Normalização pela soma dos pedaços, não por `σ_nc`.** Os intervalos particionam
`[0,1]` exatamente, então a soma dos pedaços *é* a integral — calculada com o
mesmo integrador adaptativo que calculou cada pedaço. Dividir por `σ_nc`, que num
kernel de tabela vem de um trapézio sobre os nós em `y`, misturaria dois
integradores e a identidade não fecharia. A diferença entre os dois não é
escondida: fica em `worst_norm_mismatch()`, e vale **1,04%** — é a resolução da
grade em `y`.

### `assert_composable` — o guarda oposto

`assert_comparable` exige `current` **igual**: é o guarda certo para dipolo
contra colinear.

Somar CC com NC é o caso oposto: as duas tabelas têm de ter `current` diferente,
e uma tem de ser CC e a outra NC. Somar duas tabelas de CC contaria a mesma coisa
duas vezes, e é um erro que passaria calado — o resultado seria uma seção de
choque perfeitamente plausível, duas vezes grande demais.

### Ponta a ponta

Com `σ_CC` zerada, a NC só redistribui e vaza: `Σφ + vazamento` conserva a
**4,4e-16**. Com absorção CC ligada, a regeneração aumenta o fluxo sobrevivente
em **+4,5%** sobre a mesma geometria sem regeneração.

### Gerar as tabelas

```bash
export PATH=/home/rafael/micromamba/envs/dis/bin:$PATH
export LD_LIBRARY_PATH=/home/rafael/micromamba/envs/dis/lib:$LD_LIBRARY_PATH
cd dipole
./build/sigma_nuN --model GBW --current NC --NE 300 --logEmin 3 --logEmax 14 \
    --nY 120 --use-F3 1 --pdf-set NNPDF31_nlo_as_0118      # ~95 s
```

Produz `data/sigma_nuN_NC_GBW.dat` e `data/dsigma_dy_NC_GBW.dat`, ambos com as
nove chaves e o hash do commit.
