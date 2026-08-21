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
| D5 | alto | `sigma_nuN.cpp:75` vs `scan_x.cpp`/`main.cpp` | duas convenções de α no mesmo código; o manual (§3.3) documenta a que **não** alimenta σ |
| D6 | alto | `dipole_models.hpp:6-9` + `parameters.hpp:20` | parâmetros GBW são o ajuste **sem charm** (σ₀=29,12 mb, λ=0,277, x₀=0,41e-4, que vem com m_{u,d,s}=0,14 GeV), mas o código usa m=0,03 e soma canal (c,s) |
| D7 | alto | `optical_depth.cpp:83` | `iss >> Enu >> sigma_cm2` lê a **coluna 2** (`sigma_GeV_minus2`); τ sai 2,6e27× maior |
| D8 | médio | `sigma_nuN.cpp:33` | `largeXFactor = (1-x)^7` ad hoc, aplicado a F_T/F_L/F₂ mas **não** a xF₃ |
| D9 | médio | `dipole_models.cpp:106` vs `sigma_nuN.cpp:196` | IIM zera em x ≥ 1e-2; GBW não tem corte e vai até x=0,999999 |
| D10 | médio | `parameters.hpp:87` | α do NC 2× baixo: deveria ser (g_Z/2)²/4π = √2·G_F·M_Z²/(4π) |
| D11 | médio | `parameters.hpp` | nenhum elemento de CKM; só 2 dos 6 canais W⁺→qq̄′ |
| D12 | médio | todos os `.dat` | zero proveniência — só `# model GBW` |
| D13 | médio | — | nenhuma validação contra HERA, que é o dado a que GBW e bCGC **foram ajustados** |
| D14 | baixo | `dipole_models.hpp:16` | `x0 = 0.00069e-6` do bCGC — verificar contra a fonte (padrão dos ajustes publicados é `0.00069e-4`) |
| D15 | baixo | `dipole_models.hpp:22` | `Nb = 80`, com comentário do próprio código pedindo 400 |
| D16 | baixo | `sigma_nuN.cpp:283` | `--NE 1` → divisão por `NE-1` = 0 → NaN → LHAPDF aborta |
| D17 | baixo | `integrals.cpp` | F_T e F_L em duas integrais 2D separadas, recomputando K₀/K₁ |
| D18 | baixo | `optical_depth.cpp` vs `tau_profile.cpp` | um usa `ρ/m_b`, outro `ρ·N_A` (0,7% de diferença, mas inconsistente) |
| D19 | baixo | `generate_collapsar_profile.cpp` | comentário documenta ρ=ρ_c/[1+(r/r_c)³]; código implementa lei de potência com corte exponencial. Citação [1] é uma mistura (O'Connor & Ott 2010 = CQG 27, 114103; ApJS 219, 24 = O'Connor 2015) |

---

## 4. Decisões pendentes (bloqueiam fases específicas)

### Q1 — Convenção de α  *(bloqueia F5; afeta a normalização final)*

`sigma_nuN.cpp:75-77` divide F_T e F_L por `alphaEW`. No limite EM (g_V=e_f, g_A=0, α=α_em) a
máquina reproduz o F₂^em padrão exatamente — **verificado**. Mas para CC isso deixa carga²
efetiva = g_V²+g_A² = **2** por canal, enquanto a fórmula mestre (G_F²s/2π com F₂ = 2xq) quer
peso 1×|V_qq′|².

Não afirmo que está errado: o σ convergido já sai em 0,5–0,9× Gandhi; tirar mais um fator 2 o
levaria a 0,25–0,45×, que parece baixo demais. **Precisa da tese do Alex / da referência que
fixou a convenção.**

Agravante: `scan_x.cpp` e `main.cpp` **não** dividem por α. O `data/scan_x_F2_CC_models.dat`
guarda F₂ 236× (=1/α_CC) menor que o F₂ que entra em σ, e põe essa coluna lado a lado com `xF3`
do LHAPDF na normalização padrão.

### Q2 — `x0` do bCGC  *(bloqueia a regeneração do IIM)*

`0.00069e-6` = 6,9e-10. Os ajustes Watt–Kowalski publicados usam o formato `0.00105e-4`. Um
irmão `0.00069e-4` = 6,9e-8 seria 100× maior. Como Q_s² ∝ (x₀/x)^λ, isso vale ~2,5× em σ_IIM.

### Q3 — Conjunto GBW  *(bloqueia F6)*

Escolher **um**: (a) ajuste sem charm (σ₀=29,12, λ=0,277, x₀=0,41e-4, m_q=0,14 GeV, **sem** canal
c s̄), ou (b) ajuste com charm (σ₀=23,03, λ=0,288, x₀=3,04e-4, m_{u,d,s}=0,14, m_c=1,4 GeV,
**com** canal c s̄). Hoje o código mistura os dois.

### Q4 — Política para x > 1e-2  *(bloqueia F7)*

Nem extrapolar (GBW) nem zerar (IIM). Opções: casar com PDF colinear acima de x_match, ou
declarar a tabela válida só acima de E onde x_típ < x_match e recusar abaixo.

---

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
- F1 passa a dar χ²/ponto ≲ 2.

---

### F3 — Tabular F(x,Q²) e desacoplar os laços

**Objetivo:** hoje a integral 2D roda *dentro* dos laços de x e Q² — ~9×10⁹ avaliações de Bessel.
É por isso que as grades foram apertadas até quebrar.

- Tabular F_T(x,Q²) e F_L(x,Q²) uma vez em grade log, interpolar bilinear em (ln x, ln Q²).
- Grade default: 60 pontos em x ∈ [1e-14, 1], 40 em Q² ∈ [1, 2e8] — parametrizável.
- Erro de interpolação medido e reportado, não assumido.

**Aceitação:** σ(1e6) com tabulação difere < 1% do cálculo direto convergido, e roda ≥ 50× mais
rápido. Sem isso, F4 é inviável.

---

### F4 — Grades como parâmetro de verdade

- `Makefile` passa `--Nr --Nz --NlogQ --Nlogx` explicitamente (D3).
- Remover o segundo conjunto de defaults: `integrals.hpp` e `main()` param de discordar.
- `NlogQ` adaptativo ou grade concentrada em torno de Q² ~ M_W², já que o pico não se move com E
  mas a faixa de integração cresce (D4).
- Corrigir `NE=1` (D16).

**Aceitação:** `tests/test_conditioning.py` passa — amplificação ≤ 5 (hoje: 3,55e5). σ(E)
monotônica em toda a faixa 1e3–1e14 (a tabela atual tem 29 passos decrescentes), e inclinação
log-log por década ≥ 0,25 acima de 1e7 GeV.

---

### F5 — Convenção de α  *(BLOQUEADA por Q1)*

- Fixar uma convenção, documentar em `Projeto_DIS.pdf` §3.3, aplicar aos três arquivos.
- Teste: limite EM continua exato depois da unificação.

---

### F6 — Parâmetros físicos  *(BLOQUEADA por Q3)*

- Conjunto GBW consistente (parâmetros ↔ massas ↔ canais) — D6.
- Remover ou justificar `(1-x)^7`; se ficar, aplicar também a xF₃ — D8.
- Incluir CKM e os 6 canais W⁺→qq̄′ — D11.
- Corrigir α do NC — D10.
- `Nb` do bCGC para 400 e medir o efeito — D15.

**Aceitação:** F1 (HERA) não degrada; cada mudança reportada isoladamente no seu efeito sobre σ.

---

### F7 — Política para x > 1e-2  *(BLOQUEADA por Q4)*

Mesma política para GBW e IIM, qualquer que seja. Hoje a assimetria é o que produz
GBW/IIM = 275 em 1e3 GeV e 3,3 em 1e14.

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

**Aceitação:** dentro de fator 2 do pQCD em 1e5–1e12 GeV; monotônica; proveniência completa.

---

## 6. Critério de aceitação da campanha

1. χ²/ponto ≲ 2 contra σ_red do HERA na janela de ajuste do GBW.
2. 0,05 ≤ F_L/F₂ ≤ 0,25 em x ≲ 1e-3.
3. σ(E) monotônica em 1e3–1e14 GeV, com amplificação de condicionamento ≤ 5.
4. Dobrar todas as grades muda σ em < 1%.
5. σ dentro de fator 2 do pQCD em 1e5–1e12 GeV.
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
