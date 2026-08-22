# PHASIS

## `dipole/` — seções de choque ν-N em CC pelo formalismo de dipolos de cor

Cálculo de F_T, F_L, F_2 e σ_νN(E) para corrente carregada no formalismo de
dipolos de cor, com os modelos de saturação **GBW** e **bCGC/IIM**, e sua
aplicação a profundidade óptica e raio de neutrinosfera em perfis tipo
collapsar.

Implementa Kutak & Kwieciński, *Eur. Phys. J. C* **29** (2003) 521 —
eqs. (2), (8), (9), (12), (13).

### Estado: campanha de correção em andamento

As tabelas `data/sigma_nuN_CC_{GBW,IIM}.dat` distribuídas na versão inicial
estavam **40–90× acima** de qualquer cálculo moderno de σ_νN CC. A causa foi
isolada e **não era escolha de modelo — era erro de quadratura**: integração em
`r` numa grade uniforme que não resolvia o pico em r ~ 1/ε, combinada com
quase-singularidades integráveis em z→0 e z→1 recebendo peso espúrio.

O plano completo, a evidência e o progresso estão em
**[dipole/CAMPANHA_CORRECAO.md](dipole/CAMPANHA_CORRECAO.md)**.

| Fase | Assunto | Estado |
|------|---------|--------|
| F0 | Rede de segurança: oráculo, baselines congelados | concluída |
| F1 | Validação do limite EM contra σ_red do HERA | concluída |
| F2 | Quadratura: log r + z resolvido nas duas pontas | concluída |
| F3 | Tabular F(x,Q²) e desacoplar os laços | concluída |
| F4 | Densidade de nós e suavidade de σ(E) | concluída |
| F5 | Convenção de α unificada nos três arquivos | concluída |
| F7 | Política de x grande igual para GBW e IIM | concluída |
| F8 | Camada de aplicação (índice de coluna em `optical_depth`) | concluída |
| F9/F10 | Proveniência e regeneração das tabelas | GBW concluída |
| F6 | Efeito das massas do ajuste | pendente |
| — | **IIM/bCGC falha a validação contra o HERA** | **bloqueada** |

> **A tabela `sigma_nuN_CC_IIM.dat` não deve ser usada em produção.** Com a quadratura
> consertada, o bCGC passou a dar ⟨modelo/dado⟩ = 0,10 contra σ_red do HERA — o dado a que ele
> foi ajustado — e a razão varia de 0,07 a 0,17 ao longo da janela, então não é só normalização.
> O cabeçalho da tabela traz o aviso e como reproduzir. Resolver exige Rezaeian & Schmidt,
> PRD **88** (2013) 074016. Ver seção 4, Q2 do plano.

### Compilar e testar

```bash
export CONDA_PREFIX=$HOME/micromamba/envs/dis
export PATH=$CONDA_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

cd dipole
make test-quadrature    # aceitação da F2: F_2, convergência, normalização vs KK
make test-table         # aceitação da F3/F4: interpolação, convergência, suavidade
make validate-hera      # limite EM vs 239 pontos H1+ZEUS
python3 tests/test_oracle.py
make sigma-nuN          # gera data/sigma_nuN_CC_GBW.dat
```

Dependências: g++ (C++17), boost (headers), LHAPDF 6.

### Bibliografia

Os PDFs de referência ficam em `refs/`, **fora do controle de versão** (material
sob copyright). As referências estão citadas por DOI/arXiv na seção 4 de
`dipole/CAMPANHA_CORRECAO.md`.
