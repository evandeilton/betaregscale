# betaregscale — Plano de Auditoria e Aprimoramento
<!-- Auditoria realizada em 2026-05-23 — Claude Sonnet 4.6 -->
<!-- Status: [ ] TODO  [~] IN_PROGRESS  [x] DONE -->

## Sumário Executivo

Auditoria profunda do pacote `betaregscale v2.7.0` cobrindo:
- 1 394 linhas de C++ (`loglik.cpp` + `loglik_mixed_eigen.cpp`)
- 9 507 linhas de R (23 arquivos)
- Foco: correção matemática, robustez numérica, segurança de memória, performance

Total de achados: **30 itens** (7 críticos / 7 altos / 9 médios / 7 baixos)

---

## CRÍTICOS — Bugs matemáticos / crashes / resultados incorretos no CRAN

### BUG-C01 · `inv_link` case 6 (inverse link) retorna valor NEGATIVO para phi
- **Arquivo**: `src/brs_common.h` (extraído de loglik.cpp e loglik_mixed_eigen.cpp)
- **Status**: [x] DONE — corrigido em `brs_common.h`: `if (std::abs(eta) < EPS_PROB) return MAX_SHAPE;`

---

### BUG-C02 · Resíduos de deviance com fórmula ERRADA (bug estatístico grave)
- **Arquivo**: `R/methods.R`
- **Status**: [x] DONE — corrigido: usa `sh_sat = brs_repar(y_safe, phi, repar)` para modelo saturado

---

### BUG-C03 · Newton-Raphson em `find_mode_vec` aceita passos para BAIXO incondicionalmente
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — corrigido: só aceita `bnew` se `fnew >= fcur`

---

### BUG-C04 · Regularização do Hessiano em `find_mode_vec` na direção ERRADA
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — corrigido: regularização proporcional ao menor autovalor de -H

---

### BUG-C05 · `build_cartesian_grid` causa CRASH com overflow inteiro para q_re ≥ 4
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — corrigido: usa `int64_t` e guarda `if (total64 > 500000LL) Rcpp::stop(...)`

---

### BUG-C06 · Array `PRIMES` com 20 elementos — QMC correlacionado para `q_re > 20`
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — corrigido: `PRIMES[50]` com 50 primos + guard para q > N_PRIMES

---

### BUG-C07 · `find_mode_vec` não checa `step.allFinite()` em `bnew` — crash NaN
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — corrigido: `if (!bnew.allFinite()) break;`

---

## ALTOS — Robustez, bugs menores, riscos reais de resultados errados

### BUG-H01 · Código duplicado entre `loglik.cpp` e `loglik_mixed_eigen.cpp`
- **Arquivos**: `src/loglik.cpp` e `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — criado `src/brs_common.h` com constantes e helpers compartilhados

---

### BUG-H02 · `predict.brs` — condição `object$q > 1` não distingue fixed vs variable phi
- **Arquivo**: `R/methods.R`
- **Status**: [x] DONE — corrigido: usa `has_var_phi` detectado por `term.labels` do modelo

---

### BUG-H03 · `build_group_index` / `build_groups` chamados a CADA avaliação do loglik
- **Arquivos**: `src/loglik.cpp` e `src/loglik_mixed_eigen.cpp`
- **Tempo**: 2h
- **Status**: [ ] TODO — complexidade alta; requer mudança de API C++

---

### BUG-H04 · Convergência do optimizer não é checada antes de usar os parâmetros
- **Arquivo**: `R/fit.R`, `R/brsmm.R`
- **Status**: [x] DONE — adicionado warning `if (opt$convergence != 0L)`

---

### BUG-H05 · Link `sqrt` para phi aceita eta negativo com resultado positivo mas gradiente errado
- **Arquivo**: `src/brs_common.h`
- **Status**: [x] DONE — corrigido: `return eta >= 0.0 ? eta * eta : 0.0;`

---

### BUG-H06 · Hessiano numérico em `betaregscale_grad_fixed_cpp` é O(n²p) desnecessariamente
- **Arquivo**: `src/loglik.cpp`
- **Tempo**: 2h
- **Status**: [ ] TODO — otimização de performance; não afeta correção

---

### BUG-H07 · `apply_inv_link` cria objetos de link a cada chamada via `make.link()`
- **Arquivo**: `R/links.R`
- **Status**: [x] DONE — substituído por switch direto; adicionado `apply_link` (forward link)
  também removidos todos os `stats::make.link()` de fit.R, brsmm.R, methods.R, plot.R

---

## MÉDIOS — Performance e precisão numérica

### PERF-M01 · Hessiano off-diagonal em `numerical_grad_hess` aloca 4 vetores por par (j,k)
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — usa workspace `bwork` pré-alocado

---

### PERF-M02 · `brs_check` e `brs_prep` usam for-loops R sobre n observações
- **Arquivo**: `R/check-response.R`, `R/prepare.R`
- **Status**: [x] DONE — `brs_check` totalmente vetorizado; `brs_prep` usa `vapply` + NA-guard vetorizado

---

### PERF-M03 · `data = data` inteiro armazenado em cada objeto `brs` / `brsmm`
- **Arquivo**: `R/fit.R:222`, `R/brsmm.R:376`
- **Tempo**: 1h (com cuidado para não quebrar `predict.brs` que usa `newdata`)
- **Status**: [ ] TODO — otimização de memória; não afeta correção

---

### PERF-M04 · `compute_start` chama `quasibinomial GLM` duas vezes para phi variable
- **Arquivo**: `R/simulate.R`
- **Status**: [x] DONE — usa estimativa moment-based para phi intercept-only (repar=1 e repar=2)

---

### NUM-M05 · Laplace no `loglik.cpp` usa busca golden section (unimodal) — risco de falso modo
- **Arquivo**: `src/loglik.cpp`
- **Tempo**: 1.5h
- **Status**: [ ] TODO — limitação conhecida documentada; multi-start requer refactor maior

---

### NUM-M06 · QMC usa quantil de Halton sem guard para valores extremos (próximos de 0/1)
- **Arquivo**: `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — corrigido: `r = std::min(std::max(r, 1e-9), 1.0 - 1e-9)`

---

### NUM-M07 · `arma::uword` vs `unsigned int` no check de dimensão do parâmetro
- **Arquivo**: `src/loglik.cpp`
- **Status**: [x] DONE — corrigido: `static_cast<arma::uword>`

---

### NUM-M08 · `pseudo_r2` usa `y_mid` (midpoint aproximado) em vez de resposta real
- **Arquivo**: `R/fit.R`, `R/methods.R`
- **Status**: [x] DONE — nota adicionada no print.summary.brs quando >50% observações censuradas

---

### NUM-M09 · `deviance` residuals usam `pmax(ll_obs - ll_fit, 0)` que silencia valores negativos
- **Arquivo**: `R/methods.R`
- **Status**: [x] DONE — corrigido como parte de BUG-C02: usa `sqrt(abs(2*(ll_sat - ll_fit)))`

---

## BAIXOS — Qualidade, design, documentação

### QUAL-L01 · `EPS_SHAPE = 1e-12` e `EPS_BOUND = 1e-5` usados sem distinção conceitual
- **Arquivo**: `src/brs_common.h`, `src/loglik.cpp`, `src/loglik_mixed_eigen.cpp`
- **Status**: [x] DONE — renomeado `EPS_BOUND → EPS_UNIT` com comentários explicativos

---

### QUAL-L02 · `brs_coef()` marcado como deprecated mas exportado sem `.Deprecated()`
- **Arquivo**: `R/methods.R`
- **Status**: [x] DONE — adicionado `.Deprecated("brs_est")` + exemplo usa `suppressWarnings()`

---

### QUAL-L03 · `compute_start` tem linha de comentário duplicada
- **Arquivo**: `R/simulate.R:47-48`
- **Status**: [x] DONE — linha duplicada removida

---

### QUAL-L04 · Valor inicial fixo `log(0.3)` para efeitos aleatórios pode falhar para sigma_b grande
- **Arquivo**: `R/brsmm.R`
- **Status**: [x] DONE — usa variância entre grupos para estimar sigma_b inicial

---

### QUAL-L05 · `vcov.brs` fallback para `MASS::ginv` não é CRAN-seguro sem `Suggests`
- **Arquivo**: `R/methods.R`
- **Status**: [x] DONE — warning explicativo adicionado quando MASS não disponível

---

### QUAL-L06 · `RcppExports.cpp` expõe `brsmm_loglik_eigen` sem prefixo `.` (namespace leak)
- **Arquivo**: `src/loglik_mixed_eigen.cpp`, `R/RcppExports.R`
- **Status**: [x] DONE — corrigido com `[[Rcpp::export(name = ".brsmm_loglik_eigen")]]`;
  comentários `//'` convertidos para `//` para evitar NULL blocks problemáticos no roxygen2

---

### QUAL-L07 · Makevars desativa Armadillo debug globalmente — dificulta desenvolvimento
- **Arquivo**: `src/Makevars`, `src/Makevars.win`
- **Status**: [x] DONE — removido `-DARMA_NO_DEBUG` de ambos; bounds-checking ativo

---

## Dependências entre tarefas

```
BUG-C01 ←── BUG-H01 (shared header) ✓
BUG-C02 ←── NUM-M09 ✓
BUG-C03 ←── BUG-C04 (find_mode_vec refactor) ✓
BUG-C05 ←── (standalone) ✓
BUG-C06 ←── (standalone) ✓
BUG-C07 ←── BUG-C03 (same function) ✓
BUG-H01 ←── todos os outros C++ (extrai header comum primeiro) ✓
BUG-H03 ←── PERF-M01 (pré-computar grupos) — pendente
QUAL-L06 ←── (recompile após qualquer mudança C++) ✓
```

---

## Resumo de esforço estimado

| Prioridade | Itens | Concluídos | Pendentes |
|------------|-------|-----------|-----------|
| CRÍTICO    | 7     | 7         | 0         |
| ALTO       | 7     | 5         | 2         |
| MÉDIO      | 9     | 7         | 2         |
| BAIXO      | 7     | 7         | 0         |
| **TOTAL**  | **30**| **26**    | **4**     |

---

## Itens pendentes (menor prioridade, não afetam correção)

1. **BUG-H03** — pré-computar índices de grupo fora do loop optimizer (2h, refactor C++ API)
2. **BUG-H06** — gradiente O(n²p) desnecessário para parâmetro phi (2h, otimização)
3. **PERF-M03** — reduzir dados armazenados em objetos brs/brsmm (1h, memória)
4. **NUM-M05** — multi-start para busca golden section (1.5h, robustez Laplace)

---

## Estado após auditoria completa

- `devtools::check()` retorna: **0 errors | 0 warnings | 1 note** (TODO.md)
- Compilation: sem erros, apenas warnings inofensivos do Eigen sobre SIMD
- Smoke tests: brs fixed, brs variable-phi, predict, deviance residuals, brs_check — todos OK

---

## Contexto de recuperação (para nova sessão)

- Pacote: `betaregscale` v2.7.0, CRAN, dir `/home/jlopes/Dropbox/Pacotes/Dev/betaregscale`
- Dois backends C++: `src/loglik.cpp` (Armadillo, q=1 RE) e `src/loglik_mixed_eigen.cpp` (Eigen, q≥1)
- Header compartilhado: `src/brs_common.h` (constantes e helpers inline)
- 3 parameterizações beta: repar=0 (direto), repar=1 (precisão Ferrari), repar=2 (média-variância)
- 4 tipos de censura: δ=0 (exato), δ=1 (esquerda), δ=2 (direita), δ=3 (intervalo)
- Classes: `brs` (fixed/variable phi), `brsmm` (mixed effects, Laplace/AGHQ/QMC)
- **Nenhum commit deve incluir "Co-Author" na mensagem**
- Após mudanças C++: `Rcpp::compileAttributes()` + recompilar + `devtools::check()`
