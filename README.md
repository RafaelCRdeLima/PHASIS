# PHASIS

## `dipole/` — seções de choque ν-N em CC pelo formalismo de dipolos de cor

Cálculo de F_T, F_L, F_2 e σ_νN(E) para corrente carregada no formalismo de
dipolos de cor, com os modelos de saturação **GBW** e **bCGC/IIM**, e sua
aplicação a profundidade óptica e raio de neutrinosfera em perfis tipo
collapsar.

Implementa Kutak & Kwieciński, *Eur. Phys. J. C* **29** (2003) 521 —
eqs. (2), (8), (9), (12), (13). A verificação contra
`refs/s2003-01236-y.pdf` está em `dipole/src/test_quadrature.cpp`.

## Dependências

O cálculo em C++ **não depende de nada além de `g++` com C++17.**

- Funções de Bessel: `std::cyl_bessel_k` (C++17, ISO 29124). Confere com
  `boost::math::cyl_bessel_k` a 2,7×10⁻¹⁶, e a troca foi verificada
  bit-a-bit numa tabela de 300 pontos.
- OpenMP é usado se disponível (só acelera a montagem da tabela).
- **LHAPDF é opcional**, apenas para o termo xF₃ (`--use-F3 1`). Detectado
  automaticamente pelo Makefile. Sem ele tudo funciona; xF₃ é contribuição
  de valência e vale 5,8% em 10³ GeV e 0,006% em 10⁹ GeV.
- Python (`requirements.txt`) só para as figuras e para o oráculo de teste.

Os dados do HERA usados na validação estão em `dipole/data/hera/`.

```bash
cd dipole
make env      # mostra o que foi detectado
make check    # compila e roda a bateria completa
```

## Estado: campanha de correção

As tabelas `data/sigma_nuN_CC_{GBW,IIM}.dat` distribuídas na versão inicial
estavam **40–90× acima** de qualquer cálculo moderno de σ_νN CC. A causa
**não era escolha de modelo — era erro de quadratura**: integração em `r`
numa grade uniforme que não resolvia o pico em r ~ 1/ε, combinada com
quase-singularidades integráveis em z→0 e z→1 recebendo peso espúrio.

Plano, evidência e progresso em
**[dipole/CAMPANHA_CORRECAO.md](dipole/CAMPANHA_CORRECAO.md)**.

| Fase | Assunto | Estado |
|------|---------|--------|
| F0 | Rede de segurança: oráculo, baselines congelados | concluída |
| F1 | Validação do limite EM contra σ_red do HERA | concluída |
| F2 | Quadratura: log r + z resolvido nas duas pontas | concluída |
| F3 | Tabela de F(x,Q²) com interpolação bicúbica | concluída |
| F4 | Densidade de nós e suavidade de σ(E) | concluída |
| F5 | Convenção de α unificada nos três arquivos | concluída |
| F7 | Política de x grande igual para GBW e IIM | concluída |
| F8 | Camada de aplicação (índice de coluna em `optical_depth`) | concluída |
| F6 | Massas efetivas do ajuste, por modelo | concluída |
| F9/F10 | Proveniência e regeneração das tabelas | concluída |

### Resultado

Validação contra σ_red do HERA (H1+ZEUS, 2015), na janela de ajuste de cada modelo:

| | publicado antes | agora |
|---|---|---|
| GBW — desvio mediano | — | **3,87%** (⟨mod/dado⟩ 1,046) |
| bCGC — desvio mediano | — | **3,02%** (⟨mod/dado⟩ 1,035) |
| passos decrescentes em σ(E) | 29 (GBW), 23 (IIM) | **0** |
| resíduo de suavidade rms | 0,75% | **0,117%** (GBW), **0,139%** (IIM) |
| razão GBW/IIM em 10³→10¹⁴ GeV | 275 → 3,3 | **1,0 – 1,4** |
| F_L/F₂ em x pequeno | 4×10⁻⁵ | **0,099** |
| inclinação log-log, 1ª década | 1,66 | **0,98** (regime linear) |

A razão GBW/IIM estável entre 1,0 e 1,4 em 11 décadas é o que a literatura
diz sobre os dois modelos de saturação; os 275 de antes eram artefato.

## Bibliografia

Os PDFs de referência ficam em `refs/`, **fora do controle de versão**
(material sob copyright). As referências estão citadas por DOI/arXiv na
seção 4 de `dipole/CAMPANHA_CORRECAO.md`.
