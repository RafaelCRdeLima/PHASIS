# docs/

Dois documentos, do mesmo conjunto tipográfico.

| arquivo | o que é |
|---|---|
| `relatorio.pdf` | **Equações e Pendências** — onde vive cada equação física (30, com arquivo e linha) e os 14 problemas, ordenados pelo que bloqueia a varredura de produção. |
| `comandos.pdf` | **Lista de comandos** — como rodar, e que teste cada comando permite fazer. |

Os dois se referenciam: o relatório aponta para `comandos.pdf` para
reproduzir cada número; a lista de comandos aponta para os problemas do
relatório pelo identificador (`P1`, `P4`, ...).

Recompilar:

```bash
cd docs && pdflatex relatorio.tex && pdflatex relatorio.tex
cd docs && pdflatex comandos.tex  && pdflatex comandos.tex
```

Duas passadas: a segunda resolve as referências e a numeração de página.
Não precisa de bibtex nem de pacote fora do texlive padrão.
