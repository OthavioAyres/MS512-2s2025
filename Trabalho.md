## MS512 – Trabalho: Discussão dos Resultados
**Alunos**: Othavio Henrique de Jesus Ayres RA: 246666, e  Lara Maria Herrera - RA:236181

---

### Questão 1

Notou-se que com a base original \((1, t)\), a matriz é mal condicionada, o que torna as equações normais numericamente muito instáveis, isto faz que o método QR seja mais eficiente (robusto).  Com a mudança de base e reescalonamento da matriz,temos uma matriz bem condicionada e isto faz com que o metodo QR e as equações normais cheguem no mesmo resultado. isso mostra que uma boa escolha de base pode ser o suficiente para resolver problemas numéricos.

---

### Questão 2

A questão 2 nos mostrou que os erros dos autovetores aproximados decaem geometricamente até atingir o limite da precisão de máquina (ordem de \(10^{-16}\)). Além disto, as razões de erro observadas se estabilizam em torno de \(0{,}0445\), valor que coincide muito bem com a taxa teórica \(\left|\frac{\lambda_1 - 8}{\lambda_2 - 8}\right| \approx 0{,}04453\), confirmando na prática a teoria da iteração inversa com deslocamento.

---

### Questão 3 – Comparando métodos iterativos para o problema de Poisson

#### Tabelas de resultados 

**Jacobi**

| h    | Dimensão da matriz | Iterações até convergência |
|------|--------------------|----------------------------|
| 1/10 | 81                 | 293                        |
| 1/20 | 361                | 1077                       |
| 1/40 | 1521               | 3882                       |
| 1/80 | 6241               | 13766                      |

**Gauss–Seidel padrão**

| h    | Dimensão da matriz | Iterações até convergência |
|------|--------------------|----------------------------|
| 1/10 | 81                 | 157                        |
| 1/20 | 361                | 574                        |
| 1/40 | 1521               | 2069                       |
| 1/80 | 6241               | 7364                       |

**SOR – estudo de convergência**

| ω    | 1/20 | 1/40 | 1/80  |
|------|------|------|-------|
| 0.80 | 836  | 3000 | 10638 |
| 1.00 | 574  | 2069 | 7364  |
| 1.50 | 202  | 752  | 2706  |
| 1.80 | 88   | 250  | 978   |
| 1.90 | 176  | 179  | 449   |
| 1.95 | 349  | 354  | 354   |
| 1.97 | 585  | 579  | 599   |
| 2.00 | inf  | inf  | inf   |


- O método de **Jacobi** é o mais simples, e coerentemente também é o mais lento, pois  o número de iterações cresce rapidamente conforme diminuimos 1/h. O método de **Gauss–Seidel** reduz o número de iterações em torno da metade do numero de iterações do metodo de Jacobi, oque implica em um grande salto de eficacia. Além disto, o método **SOR** (com (\omega\) entre 1,5~1,9) foi metodo mais robusto, pois reduziu significativamente o numero de itareçaões em relação a Gauss–Seidel. Entretante, com um \(\omega < 1\) este metodo piora a convergência e com \(\omega >= 2\) levou a divergência (ou a um numero de iterações extremamente alto)
