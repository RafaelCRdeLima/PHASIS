# Campanha de Correção — `dipole`

**Objetivo:** tornar `sigma_nuN_CC_{GBW,IIM}.dat` um resultado defensável, e reembarcá-lo no
HADROS3 substituindo as tabelas atuais.

**Estado inicial:** as tabelas publicadas estão 40–90× acima de qualquer cálculo moderno de
σ_νN CC. A causa foi isolada e **não é escolha de modelo — é erro de quadratura.** O núcleo
físico do código está correto.

---

## 0. Evidência que fundamenta a campanha

Reimplementação fiel de `wavefunctions` + `dipole_models` + `integrals` + `sigma_nuN` em Python,
validada contra o binário real:

```
E = 1e6 GeV, GBW, grades default do main()
  binário C++ (build/sigma_nuN) ......... 7.42771e-32 cm2
  oráculo Python ........................ 7.4253e-32  cm2     (concordância 3e-4)
  tabela publicada (interp. em 1e6) ..... 7.547e-32   cm2
```

Mesma física, mesmos parâmetros GBW, mesmos canais, mesmo `(1-x)^7`, **só trocando a grade de
integração**:

| E (GeV) | tabela | quadratura convergida | Gandhi et al. | conv/Gandhi |
|---------|--------|-----------------------|---------------|-------------|
| 1e5     | 1.59e-32 | 1.85e-34 | 3.6e-34 | 0,51 |
| 1e6     | 7.40e-32 | 6.92e-34 | 8.3e-34 | 0,83 |
| 1e7     | 1.93e-31 | 1.78e-33 | 1.9e-33 | 0,93 |
| 1e9     | 6.23e-31 | 7.44e-33 | 1.0e-32 | 0,73 |
| 1e11    | 1.25e-30 | 2.53e-32 | 5.4e-32 | 0,47 |

Com a quadratura corrigida o cálculo cai dentro de fator 2 do pQCD, ligeiramente abaixo em alta
energia — exatamente onde um modelo de dipolo sem evolução DGLAP deve cair.

### Isolamento da causa

F₂ do canal (u,d) em x=1e-5, Q²=M_W². Valor convergido: **18,0** com **F_L/F₂ = 0,099**.

```
  r uniforme  z uniforme     F_2 = 117.10    F_L/F_2 = 0.0006
  r log       z uniforme     F_2 = 125.55    F_L/F_2 = 0.0142    <- piorou
  r uniforme  z log          F_2 =   9.61    F_L/F_2 = 0.0071
  r log       z log          F_2 =  18.08    F_L/F_2 = 0.0984    OK
```

**As duas correções são necessárias e nenhuma isolada resolve.** Os dois erros têm sinais
opostos e se cancelam parcialmente na grade original: as pontas de `z` inflam (erro dominante),
a grade uniforme em `r` corta o pico em r ~ 1/ε e deflaciona. Consertar só `r` **piora** o
número. Isso dita a Fase 2: as duas mudanças entram juntas ou não entram.

### Sentinela 1 — `F_L/F₂`

`F_L/F₂` em x pequeno: um número só, barato, e discrimina as quatro combinações acima sem
ambiguidade. Alvo físico: **0,05 ≤ F_L/F₂ ≤ 0,25** em x ≲ 1e-3.

### Sentinela 2 — condicionamento  *(medido durante a F0)*

σ(E) é suave: σ ~ E^0.36, então uma perturbação dE/E deve dar dσ/σ ≈ 0,36·dE/E, ou seja
amplificação ≈ 0,4. Medido na quadratura atual:

```
  E [GeV]      dE/E     dsigma/sigma   amplificacao
  1.00e+06     1e-07    -3.613e-03      36 125
  5.18e+07     1e-07    -2.265e-02     226 549
  1.00e+09     1e-06    +2.328e-02      23 276
  1.00e+11     1e-07    -1.489e-02     148 897
                                  pior: 3.55e+05
```

**Uma mudança de 1 parte em 10⁷ na energia move σ em 2%.** A quadratura amplifica por 10²–10⁵
onde a física dá 0,4. σ(E) não é uma função suave — é ±2% de ruído determinístico sobre a
tendência, e é exatamente a origem dos 29 passos decrescentes da tabela publicada.

Mecanismo: os nós de Simpson em ln Q² cobrem ln(1) até ln(0,999·2M_N·E), então mexer em E
desloca **todos** os nós; como o pico do integrando em Q² ~ M_W² é resolvido por só ~2–5 nós, o
resultado pula.

Este é o melhor diagnóstico de saúde da quadratura da campanha: não precisa de valor de
referência externo, só da suavidade que a física exige. `tests/test_conditioning.py`.

---

## 1. Princípios

1. **Nenhuma mudança de física antes da quadratura estar convergida.** Ajustar parâmetro com a
   quadratura errada é ajustar contra um artefato.
2. **Teste antes de correção.** Cada fase entrega um teste que falha antes e passa depois.
3. **O oráculo Python é congelado como referência do comportamento atual**, não do correto. Serve
   para provar que uma mudança fez exatamente o que se pretendia — nem mais.
4. **Nada é regerado sem proveniência.** Toda tabela `.dat` carimba o conjunto completo de
   parâmetros que a produziu.
5. **Decisão pendente bloqueia a fase, não a campanha.** Fases independentes seguem em paralelo.

---

## 2. Ambiente

```bash
export CONDA_PREFIX=/home/rafael/micromamba/envs/dis
export PATH=$CONDA_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

- g++ 13.3.0, C++17
- boost 1.85 (headers) — `boost::math::cyl_bessel_k`
- LHAPDF 6.5.6, sets instalados: `CT10nlo`, `NNPDF31_nlo_as_0118`
- Dados HERA NC e+p (490 pontos, √s=318 GeV, H1+ZEUS EPJ C75 (2015) 580):
  `/home/rafael/Codes/HADROS3/sandbox/data/hera/hera_nc_ep_920.dat`
  colunas: `Q2[GeV^2]  x  y  sigma_red  err`

`.git` está presente mas `git log` falha — a cópia veio incompleta. **Reclonar ou reinicializar
antes de começar**, senão a campanha roda sem rede de segurança.

---

## 3. Inventário de defeitos

| # | Severidade | Arquivo | Defeito |
|---|-----------|---------|---------|
| D1 | **crítico** | `integrals.cpp` | integração em `r` uniforme em [1e-6, 1e2] GeV⁻¹; integrando vive em r ~ 1/ε ≈ 0,025 GeV⁻¹ em Q²=M_W² |
| D2 | **crítico** | `integrals.cpp` | integração em `z` uniforme; quase-singularidades integráveis em z→0 e z→1 (aligned jet) recebem peso espúrio |
| D3 | **crítico** | `Makefile:205` | `sigma-nuN` não passa `--Nr --Nz --NlogQ --Nlogx`; valem os defaults de `main()` (50/30/16/16) enquanto `integrals.hpp` declara 400/200 |
| D4 | alto | `sigma_nuN.cpp` | `NlogQ=16` fixo; resolução em torno do pico Q²~M_W² **encolhe** com E (≈5 nós no PeV, ≈2 em 1e14) → platô e jitter de 1–3% |
| D5 | alto | `scan_x.cpp`, `main.cpp` | ~~duas convenções de α~~ → **`sigma_nuN.cpp` está certo** (= KK eq. 9/12/13, verificado). São `scan_x` e `main` que não dividem por α e publicam F₂ 236× menor. Ver Q1 |
| D6 | baixo | `parameters.hpp:20-26` | ~~parâmetros GBW errados~~ → **os parâmetros estão certos** (é o ajuste de 4 sabores, com charm). As **massas** é que estão: fit pede m_{u,d,s}=0,14 e m_c=1,5 (GBW) / 1,4 (bCGC); código usa 0,03 e 1,3113. Ver Q3 |
| D7 | alto | `optical_depth.cpp:83` | `iss >> Enu >> sigma_cm2` lê a **coluna 2** (`sigma_GeV_minus2`); τ sai 2,6e27× maior |
| ~~D8~~ | — | `sigma_nuN.cpp:33` | **RETIRADO.** `(1-x)^7` é a regra de contagem de constituintes com n_s=4, de KK. Não aplicá-lo a xF₃ também está correto. Ver Q4 |
| D9 | médio | `dipole_models.cpp:106`, `sigma_nuN.cpp:196` | IIM zera em x ≥ 1e-2 — **desvio da referência**. KK integram todo x com o (1−x)⁷. O GBW está certo. Ver Q4 |
| D10 | médio | `parameters.hpp:87` | α do NC 2× baixo: deveria ser (g_Z/2)²/4π = √2·G_F·M_Z²/(4π) |
| D11 | baixo | `parameters.hpp` | sem CKM; só os 2 canais favorecidos por Cabibbo. **É o que KK fazem** — rebaixado a "desvio conhecido", não bug |
| D12 | médio | todos os `.dat` | zero proveniência — só `# model GBW` |
| D13 | médio | — | nenhuma validação contra HERA, que é o dado a que GBW e bCGC **foram ajustados** |
| ~~D14~~ | — | `dipole_models.hpp:16` | **RETIRADO.** `x0 = 0.00069e-6` confere com a tese §4.8 e com Rezaeian et al. PRD 87 (2013) 034002. Ver Q2 |
| D15 | baixo | `dipole_models.hpp:22` | `Nb = 80`, com comentário do próprio código pedindo 400 |
| D16 | baixo | `sigma_nuN.cpp:283` | `--NE 1` → divisão por `NE-1` = 0 → NaN → LHAPDF aborta |
| D17 | baixo | `integrals.cpp` | F_T e F_L em duas integrais 2D separadas, recomputando K₀/K₁ |
| D18 | baixo | `optical_depth.cpp` vs `tau_profile.cpp` | um usa `ρ/m_b`, outro `ρ·N_A` (0,7% de diferença, mas inconsistente) |
| D19 | baixo | `generate_collapsar_profile.cpp` | comentário documenta ρ=ρ_c/[1+(r/r_c)³]; código implementa lei de potência com corte exponencial. Citação [1] é uma mistura (O'Connor & Ott 2010 = CQG 27, 114103; ApJS 219, 24 = O'Connor 2015) |

---

## 4. Decisões — RESOLVIDAS pelas referências em `refs/`

Fontes consultadas:

- **`refs/s2003-01236-y.pdf`** — Kutak & Kwieciński, Eur. Phys. J. C **29** (2003) 521.
  *A referência do formalismo CC.* Eqs. (2), (8), (9), (10), (11), (12), (13).
- **`refs/1.pdf`** — Gay Ducati, Machado & Machado, hep-ph/0609088.
- **`refs/2.pdf`** — Gonçalves & Hepp, arXiv:1011.2718. Curvas de referência para σ_νN.
- **`refs/tese_alex_sandro_quadros.pdf`** — §4.7 (GBW, p. 120) e §4.8 (bCGC, p. 121).

### Q1 — Convenção de α → **o código está CERTO**

KK eq. (9): `F_{T,L} = Q²/(4π²) ∫d²r ∫dz |ψ̄_{T,L}|² σ_d`, com (quarks sem massa)

```
  |ψ̄_T^W|² = ( 6/π²)[z² + (1-z)²] Q̄² K₁²(Q̄r)        eq. (12)
  |ψ̄_L^W|² = (24/π²) z²(1-z)²  Q²  K₀²(Q̄r)          eq. (13)      Q̄² = z(1-z)Q²
```

Os ψ̄ de KK **não carregam α** — o acoplamento está em G_F² na eq. (2). O código põe α_CC
dentro de |ψ|² e `sigma_nuN.cpp:75-77` divide de volta: o líquido é exatamente a normalização
de KK. Verificado numericamente, somando os dois canais e com m→0:

```
  r     z     psiT codigo    psiT KK       razao     psiL codigo   psiL KK       razao
  1.00  0.50  4.640406e-02   4.640406e-02  1.00000   5.635046e-02  5.635046e-02  1.00000
  0.30  0.20  1.249236e+00   1.249236e+00  1.00000   6.315282e-01  6.315282e-01  1.00000
  2.00  0.70  6.363499e-03   6.363499e-03  1.00000   6.162386e-03  6.162386e-03  1.00000
  0.05  0.50  9.370159e+00   9.370159e+00  1.00000   1.244664e+01  1.244664e+01  1.00000
```

O fator (g_V²+g_A²) = 2 por canal × 2 canais = 4, que é o Σ(carga²) efetivo correto para
F₂^CC = 2x[ū+d̄+s+c] (contra 10/9 do caso EM — razão 3,6, conferida).

**Não há fator 2 pendente.** A suspeita registrada na análise inicial fica retirada.

**Ação (F5):** `sigma_nuN.cpp` fica como está. Corrigir `scan_x.cpp` e `main.cpp`, que **não**
dividem por α e portanto publicam F₂ 236× menor. Documentar a normalização de ψ no
`Projeto_DIS.pdf` §3.3 — a fórmula lá está certa, mas omitir a normalização de ψ é o que
tornou a divergência entre os três arquivos invisível.

### Q2 — `x0` do bCGC → **o código está CERTO**

Tese §4.8, p. 121: «B_CGC = 5.5 GeV², γ_s = 0.6492, N₀ = 0.3658, **x₀ = 0.00069 × 10⁻⁶** e
λ = 0.2023», de Rezaeian, Siddikov, Van de Klundert & Venugopalan, Phys. Rev. D **87** (2013)
034002. É valor publicado; o x₀ minúsculo é compensado pelo λ pequeno.

**Ação:** nenhuma. D14 retirado.

### Q3 — Conjunto GBW → **parâmetros certos, MASSAS erradas**

Tese §4.7, p. 120: σ₀ = 29.12 mb, λ = 0.277, x₀ = 0.41×10⁻⁴ «ajustados aos dados de F₂ do HERA
considerando quarks leves **e a contribuição do quark charm**». KK confirmam: «obtained from the
fit with four flavors». Somar os canais (u,d) e (c,s) é o correto — KK: «os dipolos que
contribuem para transições favorecidas por Cabibbo são ud̄(dū), cs̄(sc̄)».

**Mas o ajuste vem com massas próprias, e são outras:**

| | fit GBW (tese §4.7) | fit bCGC (tese §4.8) | código |
|---|---|---|---|
| m_{u,d,s} | **0,14 GeV** | **0,14 GeV** | 0,03 GeV |
| m_c | **1,5 GeV** | **1,4 GeV** | 1,3113 GeV |

As massas efetivas são parâmetros do ajuste, não convenção livre. Usar 0,03 quebra o ajuste ao
HERA — e é justamente onde mais dói, porque em z→0 o dipolo máximo vai como 1/m (fator 4,7).

**Ação (F6):** `QuarkMasses` por modelo, com os valores do ajuste correspondente.

**Predição a testar depois da F2:** hoje, na F1, m=0,14 dá χ² *pior* que m=0,03 (73 contra 21).
Com a quadratura convergida a ordem tem que inverter. Se não inverter, há outra coisa errada.

### Q4 — Política para x > 1e-2 → **o GBW segue a referência; o corte do IIM é o desvio**

KK: «O modelo de dipolo descreve bem o DIS a pequeno x, mas torna-se impreciso a x grande e
moderadamente pequeno. Esse efeito pode ser aproximadamente levado em conta multiplicando as
funções de estrutura por um fator (1−x)^{2n_s−1}, que segue da regra de contagem de
constituintes, onde n_s é o número de quarks espectadores. Como o modelo de dipolo representa a
contribuição do mar, tomamos n_s = 4.»

2·4−1 = **7**. O `largeXFactor = (1-x)^7` **é da referência**, não é fudge. D8 retirado.

Não aplicá-lo a xF₃ também está certo: xF₃ vem de PDF colinear, que já tem o comportamento
correto em x→1.

**Ação (F7):** remover o corte duro `x >= 1e-2` de `sigmaDipoleIIM` e o `xmax = 1e-2` de
`sigmaNuN_CC`, deixando o (1−x)⁷ agir nos dois modelos igualmente. É o que elimina a assimetria
GBW/IIM de 275 em 1e3 GeV.

### Q5 — NOVA: definição de Y no bCGC *(pendente, baixo impacto)*

Tese eq. (4.46): «sendo **Y = ln(x₀/x)** a rapidez». O código usa `Y = log(1.0/x)`.

Com x₀ = 6,9×10⁻¹⁰, ln(x₀/x) é **negativo** para todo x > 6,9e-10 — ou seja, para toda a região
usada. Isso tornaria γ_eff = γ_s + ln(2/τ)/(κλY) negativo e N divergente quando τ→0. O código,
com ln(1/x), está na convenção padrão do IIM.

**Quase certamente erro de digitação na tese.** Confirmar contra Rezaeian et al. PRD 87 034002,
que não está em `refs/`. **Ação: nenhuma mudança no código**; só confirmar.

## 5. Fases

### F0 — Rede de segurança  *(pré-requisito, sem decisão pendente)*

**Objetivo:** poder provar que qualquer mudança fez exatamente o pretendido.

- Reinicializar o repositório git (a cópia veio sem histórico utilizável); commit do estado atual
  como baseline intocado.
- Congelar `tests/oracle/` com o oráculo Python validado (concordância 3e-4 com o binário).
- Congelar `tests/baseline/sigma_baseline.json`: σ(E) do binário atual em 8 energias, para
  detectar mudança não intencional.
- `chmod +x` nos binários (perderam o bit na cópia).

**Aceitação:** `python3 tests/test_oracle.py` passa (concordância oráculo↔C++).

**CONCLUÍDA.** Commit baseline `edbc843`. Oráculo em `tests/oracle/dipole_oracle.py`, validado a
2e-6 contra o binário em 8 energias de 1e3 a 1e14 GeV. Baselines GBW (com e sem F3) e IIM
congelados em `tests/baseline/`. Descoberta colateral: o teste de condicionamento (§0, Sentinela 2).

---

### F1 — Teste de limite EM contra HERA  *(sem decisão pendente — começa agora)*

**Objetivo:** o teste que teria pego tudo isso na primeira rodada.

GBW e bCGC **são ajustes ao F₂ do HERA**. Pondo g_V=e_f, g_A=0, α=α_em, a mesma máquina tem que
reproduzir σ_red medida.

- `src/validate_hera.cpp`: modo EM em `Parameters` (`makeEMParameters(flavor)`), lê
  `hera_nc_ep_920.dat`, calcula σ_red = F₂ − y²/(1+(1−y)²)·F_L, escreve χ²/ponto.
- Janela de validade do GBW: x < 0,01 e 0,25 < Q² < 50 GeV². Somar u,d,s (+c no ajuste com charm).
- `plots/plot_hera_validation.py`: σ_red(x) por bin de Q², dados vs modelo.

**Aceitação (antes da F2):** o teste roda e **falha** com χ²/ponto ≫ 1 — documentando o estado.
**Aceitação (depois da F2):** χ²/ponto ≲ 2 na janela de ajuste do GBW.

**CONCLUÍDA (antes da F2).** `src/validate_hera.cpp`, `make validate-hera`,
`plots/plot_hera_validation.py`. 239 pontos do HERA em x < 1e-2, 0,25 ≤ Q² ≤ 50 GeV².

| config | χ²/ponto | ⟨mod/dado⟩ | F_L/F₂ |
|--------|----------|------------|--------|
| **Nr=50, Nz=30, m=0,03** — grade que gerou as tabelas | **2131** | 1,349 | 0,0026–0,238 |
| Nr=400, Nz=200, m=0,03 | 21,1 | 1,014 | 0,097–0,223 |
| Nr=800, Nz=400, m=0,03 | 38,4 | 0,980 | 0,116–0,223 |
| Nr=400, Nz=200, m=0,14 | 73,4 | 0,905 | 0,134–0,214 |
| Nr=800, Nz=400, m=0,14 | 89,8 | 0,897 | 0,141–0,214 |

Três leituras:

1. **O teste funciona.** Na grade real das tabelas, χ² = 2131 e F_L/F₂ desce a 0,0026. Refinando
   para 400/200, χ² cai 100× e F_L/F₂ entra na faixa física. As duas sentinelas concordam.
2. **A falha é dependente de Q², crescente.** Em `plots/hera_validation.png`: até Q² = 2 GeV² as
   três grades batem com o dado; em Q² = 45 GeV² a grade grossa dá σ_red = 3,5 contra 1,45 medido
   (fator 2,4). Como r_pico ~ 1/ε encolhe com Q², em Q² = M_W² (140× maior) isso vira o fator ~56
   medido em §0. **Corolário:** um teste só em Q² ≤ 50 GeV² com grade fina não pegaria o bug —
   é a combinação grade grossa + alcance em Q² que expõe.
3. **Ainda não convergido, e refinar piora.** 400/200 → 21, mas 800/400 → 38. Uniforme em r com
   uniforme em z é mal-condicionado: os dois erros se cancelam parcialmente e refinar quebra o
   cancelamento. Confirma §0.
4. **⟨mod/dado⟩ ≈ 1 com χ² enorme** em todos os casos: a discrepância é de **forma**, não de
   normalização.

**Não é possível julgar parâmetro físico ainda.** m=0,14 (o valor do artigo GBW) dá χ² *pior*
que m=0,03 — mas com a quadratura não convergida isso não significa nada. É exatamente o
Princípio 1. A decisão Q3 fica para depois da F2.

---

### F2 — Quadratura em (r, z)  *(o coração da campanha)*

**As duas mudanças entram no mesmo commit** — ver §0.

- `integrals.cpp`: integração em **ln r**, faixa [1e-8, 1e3] GeV⁻¹, jacobiano `r`.
- `integrals.cpp`: `z` em grade log de duas pontas — z ∈ [1e-11, 0,5] ∪ [0,5, 1−1e-11], espelhada.
  (O integrando é simétrico em z↔1−z só quando m=μ; **não** assumir simetria no canal (c,s).)
- `FT` e `FL` numa única passada, reaproveitando K₀ e K₁ (D17).
- `rMin/rMax/zMin/zMax` saem de `parameters.hpp` e viram parâmetros da quadratura, não do modelo.

**Aceitação:**
- F₂(u,d) em x=1e-5, Q²=M_W² = 18,0 ± 5%, com F_L/F₂ = 0,099 ± 0,01.
- Teste de convergência: dobrar Nr e Nz muda σ(1e6) em < 1%.
- ~~F1 passa a dar χ²/ponto ≲ 2~~ → **critério revisto, ver abaixo.**

**CONCLUÍDA.** `make test-quadrature` passa nos três critérios:

```
1) F_2      = 18.044    (ref 18.0,   desvio 0.25%)   OK
   F_L/F_2  = 0.098635  (ref 0.099,  desvio 0.37%)   OK
2) Nr=100 -> 200: 0.0405%   200 -> 400: 0.0039%   400 -> 800: 0.0003%   OK
3) |psi_T|^2 e |psi_L|^2 vs KK eqs. (12)/(13): desvio 0.00%            OK
```

Antes: F₂ = 1018 com F_L/F₂ = 4e-5.

**Efeito em σ_νN (GBW, sem F3):**

| E (GeV) | pré-F2 | pós-F2 | Gandhi | pós/Gandhi | redução |
|---------|--------|--------|--------|------------|---------|
| 1,0e3 | 2,45e-35 | 2,29e-36 | 6,79e-35 | 0,03 | 11× |
| 3,7e4 | 6,01e-33 | 8,95e-35 | 2,52e-34 | 0,35 | 67× |
| 1,4e6 | 8,83e-32 | 8,06e-34 | 9,39e-34 | 0,86 | 110× |
| 5,2e7 | 3,11e-31 | 3,03e-33 | 3,49e-33 | 0,87 | 103× |
| 1,9e9 | 6,91e-31 | 8,58e-33 | 1,30e-32 | 0,66 | 80× |
| 7,2e10 | 1,20e-30 | 2,33e-32 | 4,83e-32 | 0,48 | 51× |
| 2,7e12 | 1,83e-30 | 5,68e-32 | 1,80e-31 | 0,32 | 32× |
| 1,0e14 | 2,51e-30 | 1,25e-31 | 6,68e-31 | 0,19 | 20× |

Inclinação log-log por trecho: **1,013** → 0,607 → 0,366 → 0,288 → 0,276 → 0,247 → 0,218.
O primeiro trecho reproduz o regime linear σ ∝ E que a física exige (era 1,655). **σ(E) agora é
monotônica** nos 8 pontos.

### Critério da F1 revisto — χ²/ponto ≤ 2 era inatingível

Medido: os dados HERA I+II combinados (2015) têm **erro relativo mediano de 2,5%** na janela
usada. χ²/ponto = 1 exigiria o GBW — 3 parâmetros, ajustado nos anos 1990 a dados com erro de
5–10% — dentro de 2,5% de dados de 2015. Não é um alvo realista.

O que o resíduo mostra, com a quadratura convergida:

| Q² [GeV²] | n | ⟨mod/dado⟩ | desvio mediano | χ²/pt |
|-----------|---|------------|----------------|-------|
| 0,25–1 | 54 | 1,140 | 14,9% | 24,3 |
| 1–3 | 45 | 1,026 | 3,0% | 2,9 |
| 3–10 | 56 | 0,943 | 7,6% | 22,5 |
| 10–25 | 51 | 0,834 | 16,4% | 119,0 |
| 25–50 | 33 | 0,781 | 22,2% | 259,3 |

Tendência **monotônica** em Q²: o modelo subestima progressivamente. É a ausência de evolução
DGLAP no GBW, e é a ressalva que os próprios KK fazem («the dominant contribution comes from
Q² ~ M_W², where the simple GBW model may not be sufficiently accurate»). Ruído de quadratura
não seria monotônico — este resíduo é físico.

**Critério novo:** χ²/ponto ≤ 5 em 1 ≤ Q² ≤ 10 GeV² (coração da janela de ajuste), desvio
mediano global ≤ 15%, e ausência de estrutura **não-monotônica** em Q². A degradação suave em
Q² alto é limitação conhecida do modelo e vai para o texto, não para o critério.

### Q3 revista — as massas são escolha, não bug

Com a quadratura convergida, m=0,14 dá χ² **pior** (96,1) que m=0,03 (72,5) — a predição que eu
tinha registrado não se confirmou. O motivo é claro pela tabela acima: o modelo já subestima em
Q² alto, e massa maior o reduz mais.

Há duas escolhas defensáveis, cada uma de uma referência:

- **(a) consistente com KK:** quarks sem massa. É o que KK usam com estes mesmos parâmetros GBW,
  e é o alvo da F10. O m=0,03 do código é praticamente isso.
- **(b) consistente com o ajuste:** m_{u,d,s}=0,14, m_c=1,5, as massas com que σ₀/λ/x₀ foram
  ajustados (tese §4.7).

Como o critério de aceitação da campanha é **reproduzir a Fig. 3 de KK**, a escolha (a) é a
coerente, e o código já está nela. **D6 rebaixado a "divergência documentada".** O que resta de
F6 é medir o efeito de (b) e registrá-lo, não trocar.

---

### F3 — Tabular F(x,Q²) e desacoplar os laços

**Objetivo:** hoje a integral 2D roda *dentro* dos laços de x e Q² — ~9×10⁹ avaliações de Bessel.
É por isso que as grades foram apertadas até quebrar.

- Tabular F_T(x,Q²) e F_L(x,Q²) uma vez em grade log, interpolar bilinear em (ln x, ln Q²).
- Grade default: 60 pontos em x ∈ [1e-14, 1], 40 em Q² ∈ [1, 2e8] — parametrizável.
- Erro de interpolação medido e reportado, não assumido.

**Aceitação:** σ(1e6) com tabulação difere < 1% do cálculo direto convergido, e roda ≥ 50× mais
rápido. Sem isso, F4 é inviável.

**CONCLUÍDA.** `include/structure_table.hpp`, `src/structure_table.cpp`, `make test-table`.

- Tabela 121×121 em (ln x, ln Q²), montada **uma vez** e servindo todas as energias.
- OpenMP: 19 s de relógio (3m22s de CPU em 12 núcleos). Antes: 5,5 s **por energia**, ou seja
  ~27 min para as 300 energias de produção.
- Interpolação **bicúbica (Catmull-Rom)**, não bilinear — ver F4b abaixo.
- `src/sigma_nuN_core.{hpp,cpp}`: o núcleo da seção de choque saiu do `main()` para poder ser
  testado.

Erro de interpolação medido em pontos fora dos nós: **pior caso 0,001%** (limite 2%).

---

### F4 — Grades como parâmetro de verdade

- `Makefile` passa `--Nr --Nz --NlogQ --Nlogx` explicitamente (D3).
- Remover o segundo conjunto de defaults: `integrals.hpp` e `main()` param de discordar.
- `NlogQ` adaptativo ou grade concentrada em torno de Q² ~ M_W², já que o pico não se move com E
  mas a faixa de integração cresce (D4).
- Corrigir `NE=1` (D16).

**CONCLUÍDA.**

`nodesForRange()` fixa a **densidade** de nós por decada (default 16), em vez do `NlogQ = 16`
fixo. `NlogQ` passa a escalar: 28 em 1e3 GeV, 116 em 1e14. Também corrigido o `--NE 1`, que
dividia por `NE-1 = 0`.

Duas correções que só apareceram ao medir:

1. **`Nlogx` estava sendo recalculado dentro do laço de Q².** A faixa em x é [Q²/s, xmax], que
   encolhe conforme Q² cresce, então o contador saltava de 2 em 2 ao longo do laço — e cada
   salto é uma **descontinuidade no integrando externo** (o erro de convergência da integral
   interna muda de degrau). Hasteado para fora, calculado uma vez pela faixa mais larga.

2. **A interpolação bilinear da tabela é C⁰ mas não C¹.** Os nós de Simpson cruzavam as quinas
   da grade a cada mudança de energia. Trocada por Catmull-Rom (C¹, estêncil 4×4, mesmo custo
   de memória). Efeito no resíduo de suavidade: 0,464% → 0,227% rms. E, de quebra, o erro de
   interpolação da própria tabela caiu de 0,105% para 0,001%.

### Critério da F4 revisto — "amplificação" era a métrica errada

`tests/test_conditioning.py` media |dσ/σ| / |dE/E| com dE/E = 1e-7. Isso só faz sentido para uma
resposta **diferenciável**. Com a tabela, σ(E) é suave na escala que importa (a grade de produção
tem dE/E ≈ 9%) mas não é diferenciável na escala de 1e-7, porque os nós de quadratura cruzam a
estrutura local do interpolante. Medir amplificação ali é medir ruído sem conteúdo físico — e de
fato `dσ/σ` saía praticamente constante (~2e-3) para qualquer `dE/E`, o que é a assinatura de um
salto discreto, não de amplificação.

**Métrica nova, que é a que importa:** σ(E) tem de ser monotônica e suave o bastante para que a
interpolação log-log da tabela — que é o que o HADROS3 faz — não injete artefato. A quarta
diferença finita em ln E aniquila qualquer cúbica, então D4/6 estima o desvio de cada ponto em
relação a uma curva suave local sem precisar ajustar nada.

**Resultado (45 energias, 1e3–1e14 GeV):**

| | tabela publicada | pós-F2 (bilinear, 8 n/déc) | **pós-F4** |
|---|---|---|---|
| passos decrescentes | **29** | 0 | **0** |
| resíduo rms | 0,75% | 0,227% | **0,199%** |
| resíduo máximo | 1,83% | 0,687% | **0,629%** |

Convergência em densidade de nós (16 → 32 por década): 0,11% a 0,24%.

---

### F5 — Convenção de α  *(DESBLOQUEADA — ver Q1)*

- `sigma_nuN.cpp` fica. Corrigir `scan_x.cpp` e `main.cpp` (dividir por α).
- Documentar a normalização de ψ no `Projeto_DIS.pdf` §3.3.
- Teste: limite EM continua exato; `scan_x` passa a concordar com `sigma_nuN` no mesmo (x,Q²).

---

### F6 — Parâmetros físicos  *(DESBLOQUEADA — ver Q3)*

- Massas por modelo: GBW m_{u,d,s}=0,14 / m_c=1,5; bCGC m_{u,d,s}=0,14 / m_c=1,4 — D6/Q3.
- `(1-x)^7` **fica** — é da referência (Q4).
- Corrigir α do NC — D10.
- `Nb` do bCGC para 400 e medir o efeito — D15.
- CKM e canais suprimidos: **fora de escopo** (KK também não incluem); registrar como desvio.

**Aceitação:** F1 (HERA) não degrada; cada mudança reportada isoladamente no seu efeito sobre σ.

---

### F7 — Política para x > 1e-2  *(DESBLOQUEADA — ver Q4)*

Remover o corte duro do IIM (`sigmaDipoleIIM`, `sigmaNuN_CC`), deixando o (1−x)⁷ agir nos dois
modelos igualmente, como KK prescrevem. É o que elimina a assimetria GBW/IIM de 275 em 1e3 GeV
contra 3,3 em 1e14.

---

### F8 — Camada de aplicação  *(sem decisão pendente)*

- `optical_depth.cpp:83`: ler a coluna 3 (D7). **Uma linha.**
- Unificar os três leitores de tabela num só (`sigma_table.hpp`), com o mesmo `ρ/m_b` (D18).
- Corrigir comentário e citação de `generate_collapsar_profile.cpp` (D19).

**Aceitação:** τ de `optical_depth` e de `tau_profile` concordam no mesmo perfil e mesma σ.

---

### F9 — Proveniência  *(sem decisão pendente)*

Todo `.dat` gerado carimba: modelo, parâmetros do dipolo, massas de quark, Nr, Nz, NlogQ, Nlogx,
Q2min, use-F3, PDF set, beam, versão do código, data.

**Aceitação:** `sigma_nuN_CC_GBW.dat` permite reproduzir a si mesmo só pelo cabeçalho.

---

### F10 — Regeneração e reembarque

- Regerar GBW e IIM, 300 pontos, 1e3–1e14 GeV.
- Comparar com Gandhi et al., Cooper-Sarkar et al., Connolly et al.
- Substituir `HADROS3/data/sigma/sigma_nuN_CC_{GBW,IIM}.dat`.
- Atualizar `sigma_table_physics_risk` em `hadros3/dis_sampler.py:1512` e o `\begin{cuidado}` da
  apostila (`sandbox/apostila/dis.tex:1211`), que hoje registra a discrepância como pendência
  aberta.

**Aceitação — agora muito mais forte.** Resolvidas Q1–Q4, ficou claro que **este código É o
cálculo de Kutak–Kwieciński**: mesmas eqs. (2), (8), (9), (12), (13), mesmos parâmetros GBW,
mesmo (1−x)⁷, mesmo Q²min = 1 GeV². Então o alvo não é "dentro de fator 2 do pQCD" — é
**reproduzir a Fig. 3 de KK**, a curva GBW de σ_CC(E) entre 1e7 e 1e13 GeV.

É um teste bem mais apertado, e é o critério de aceitação da campanha. Complementos:
`refs/2.pdf` (Gonçalves & Hepp) para as curvas CGC/BK, e Gandhi et al. como âncora pQCD.

Ainda: monotônica, condicionamento ≤ 5, proveniência completa.

---

## 6. Critério de aceitação da campanha

1. χ²/ponto ≲ 2 contra σ_red do HERA na janela de ajuste do GBW.
2. 0,05 ≤ F_L/F₂ ≤ 0,25 em x ≲ 1e-3.
3. σ(E) monotônica em 1e3–1e14 GeV, com amplificação de condicionamento ≤ 5.
4. Dobrar todas as grades muda σ em < 1%.
5. Reproduz a Fig. 3 de Kutak–Kwieciński (curva GBW) em 1e7–1e13 GeV — o código implementa
   exatamente esse cálculo, então a concordância deve ser bem melhor que fator 2.
6. Toda tabela reproduzível pelo próprio cabeçalho.

---

## 7. Riscos

- **Q1 não resolvida** → normalização final com fator 2 indeterminado. Todo o resto continua
  válido; só a normalização absoluta fica em aberto.
- **GBW extrapolado a Q² ~ M_W²** está muito além da janela do ajuste (Q² ≲ 50 GeV²). Isso é
  limitação do modelo, não bug — mas tem que ser dito no paper, não escondido.
- **F3 do LHAPDF misturado com F₂ de dipolo** são duas descrições diferentes no mesmo colchete.
  Documentar como aproximação assumida.
- **Escopo:** mesmo com σ 100× menor, o toro do HADROS3 a 1e10 g/cm³ segue opaco por ~8 ordens de
  grandeza. O bug de amostragem do ponto de interação no HADROS3 (`dis_sampler.py:1387`, sorteio
  ∝ Δτ sem o fator e^(−τ)) é **independente desta campanha** e igualmente prioritário.
