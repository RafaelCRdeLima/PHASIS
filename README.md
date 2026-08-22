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
| F9/F10 | Proveniência e regeneração das tabelas | GBW concluída |
| F6 | Efeito das massas do ajuste | pendente |
| — | **IIM/bCGC falha a validação contra o HERA** | **bloqueada** |

### Resultado da tabela GBW

| | publicada | agora |
|---|---|---|
| passos decrescentes em σ(E) | 29 | **0** |
| resíduo de suavidade rms | 0,75% | **0,117%** |
| σ/Gandhi em 10⁷ GeV | ~100× | **0,95** |
| inclinação log-log, 1ª década | 1,66 | **0,98** (regime linear) |
| F_L/F₂ em x pequeno | 4×10⁻⁵ | **0,099** |

### Aviso sobre a tabela IIM

**`sigma_nuN_CC_IIM.dat` não deve ser usada em produção.** Com a quadratura
consertada, o bCGC passou a dar ⟨modelo/dado⟩ = 0,10 contra σ_red do HERA — o
dado a que ele foi ajustado — e a razão varia de 0,07 a 0,17 ao longo da
janela, então não é só normalização. O cabeçalho da tabela traz o aviso e como
reproduzir. Resolver exige Rezaeian & Schmidt, PRD **88** (2013) 074016.
Ver seção 4, Q2 do plano.

## Bibliografia

Os PDFs de referência ficam em `refs/`, **fora do controle de versão**
(material sob copyright). As referências estão citadas por DOI/arXiv na
seção 4 de `dipole/CAMPANHA_CORRECAO.md`.
