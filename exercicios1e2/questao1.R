#This exercise demonstrates that the QR method is superior to the normal
#equations method when the condition number is bad. It also shows how a change of
#basis can improve the condition number. Consider the following data.

#Conjunto de dados:
set.seed(12345)
options(digits = 15)
t <- c(0.98765431, 0.98765432, 0.98765433, 0.98765434, 0.98765435, 0.98765436, 0.98765437, 0.98765438, 0.98765439)
y <- c(rep(2.1, 3), 2, rep(1.9,2), 1.8, rep(1.7,2))
data <- data.frame(t,y)
data
#------------------------------------------------------------------------------------------
#a)Set up an over-determined system Ax = b for a linear polynomial. Use the
#basis phi (t) = 1, phi2 (i) = t. Using MATLAB, find the condition number of A.
#Find the least-squares solution by the QR decomposition method (e.g. x=A\b in MATLAB). 
#Calculate the norm of the residual.
#------------------------------------------------------------------------------------------
#Definindo o sistema Ax = b: phi1(t) = 1 e phi2(t) = t.
# Para facilitar a nomenclatura, vamos usar o sistema "X beta = Y", em que:
# X: é a matriz A 
# beta: é a solução do sistema (x)
# Y: correspondente a b

#Definindo a matriz X
X <- cbind(const = 1, t = t)
X # Os numeros de t sao tao pequenos e proximos, que praticamente a coluna 1 é uma combinacao linear da 2. Teremos problemas!

#verificando o determinante:
XTX <- t(X) %*% X
det(XTX) #o determinante é extremamente pequeno. É quase zero!  (lembrete: Determinante = 0 não tem inversa). 

#Número de condição de X. Lembrete: Kappa(A) = ||A||||A^-1||
KX <- kappa(X)
KX # A matriz é mal condicionada. o kappa(A) é gigante!

#Resolvendo agora os mínimos quadrados usando a fatoração QR: 
#A função lm() usar QR por trás. 
#Com lm()
solucao1 <- lm(y ~ X[, 2])
coef(solucao1)
#QR diretamente:
QR <- qr(X)
solucao2 <- qr.coef(QR, y)
solucao2
# NAO DEU CERTO!!!!!!!!!!!!!!!!!!!!!!!!!!!!




# Resolvendo: vamos reescalar os valores da X:
t_mean <- mean(t)
t_scaled <- (t - t_mean) * 1e6   # escalonamos e montamos nova matriz A
X_nova <- cbind(const = 1, t = t_scaled)

#numero de condição da matriz nova
kappa(X_nova) #melhorou muito

#nova decomposição QR:
QR2 <- qr(X_nova)
#pegando a Q:
Q <- qr.Q(QR2)
Q
#pegando a R:
R<-qr.R(QR2)
R

# resolver:   X = Q R 
#             X beta = y
#           (Qtransp X) beta = Qtransp y
#             R beta = Qtransp y
solucao.betas <- backsolve(R, t(Q) %*% y)   
solucao.betas #sao do sistema transformado (scaled)

beta0_scaled <- solucao.betas[1]; beta1_scaled <- solucao.betas[2]
# converter para a escala original:
# temos: y = beta0_scaled + beta1_scaled * [t_scaled]
#        y = beta0_scaled + beta1_scaled * [(t - t_mean) * 1e6]
#        y = beta0_scaled + beta1_scaled * [t * 1e6 - t_mean * 1e6]
#        y = beta0_scaled + beta1_scaled * (t * 1e6) - beta1_scaled (t_mean * 1e6)
# organizando:
#        y = (beta0_scaled - beta1_scaled * t_mean * 1e6) + (beta1_scaled * 1e6) * t

beta1_original <- beta1_scaled * 1e6
beta0_original <- beta0_scaled - (beta1_original * t_mean)

# Solução do sistema:
betas_qr <- c(intercept = beta0_original, slope = beta1_original)
betas_qr #sao do sistema original

#Resíduo e sua norma 2:
r <- y - X %*% betas_qr
r
res_norm <- norm(r, "2")
res_norm

#------------------------------------------------------------------------------------------
#b)Obtain the normal equations ATAx = ATb. Find the condition number of
#the coefficient matrix. Solve the normal equations. Calculate the norm of the
#residual (|| b - Ax||2).
#------------------------------------------------------------------------------------------
# ATAx = ATb, ou seja, na nossa nomenclatura: Xtransp X beta = Xtransp y
X
XTX <- t(X) %*% X
XTX
XTY <- t(X) %*% y
XTY

# número de condição de XTX
kappa(X)
kappa(XTX) #Muito grande, maior que a de X que já era mal condicionada. Ou seja, XTX é mal condicionada

# solução via equações normais
# resolver:   Xtransp X beta = Xtransp y 
beta_eqnormais <- solve(XTX, XTY)
beta_eqnormais

# Norma 2 do resíduo
res_eqnormais <- y - X %*% beta_eqnormais
res_eqnormais_norm <- norm(res_eqnormais, type = "2")
res_eqnormais_norm
#------------------------------------------------------------------------------------------
#c)Compare your solutions from parts (a) and (b). Plot the two polynomials and
#the data points on a single graph.
#------------------------------------------------------------------------------------------
betas_qr
beta_eqnormais

comparacao <- data.frame(
  metodo = c("QR", "NormalEQ"),
  beta0 = c(betas_qr["intercept"], beta_eqnormais["const", 1]),
  beta1 = c(betas_qr["slope"], beta_eqnormais["t", 1]),
  row.names = NULL)
comparacao



beta0_qr <- comparacao$beta0[1]
beta1_qr <- comparacao$beta1[1]

beta0_ne <- comparacao$beta0[2]
beta1_ne <- comparacao$beta1[2]

library(ggplot2)
library(dplyr)

data

df_lines <- data.frame(
  t = t,
  y = y,
  y_pred = c(beta0_qr + beta1_qr * t, beta0_ne + beta1_ne * t),
  metodo = rep(c("QR", "Normal Equations"), each = length(t))
)
df_lines

grafico1 <- df_lines %>% ggplot() +
  geom_point(aes(x = t, y = y), color = "black") +
  geom_line(aes(x = t, y = y_pred, color = metodo, linetype = metodo), size = 1.2) +
  scale_color_manual(values = c("QR" = "blue", "Normal Equations" = "red")) +
  scale_linetype_manual(values = c("QR" = "solid", "Normal Equations" = "dashed")) +
  labs(
    title = "Comparação: ajuste Eq lineares  QR vs Equações Normais",
    x = "t",
    y = "y"
  ) +
  theme_minimal(base_size = 14)
grafico1
#Pelo gráfico, vemos que a reta via cecomposição QR é mais estável.

#------------------------------------------------------------------------------------------
#d)Now set up a different over-determined system Âˆx = ˆb for the same data using
#the basis phi (t) = 1, phi2 (i) = 3*10^7*(t-0.98765435). Find the condition
#number of A. Find the least-squares solution by the QR decomposition method.
#Calculate the norm of the residual.
#------------------------------------------------------------------------------------------

# Vamos estabelecer a nova matriz X 
t_novo <- 3*(10^7)*(t-0.98765435)
t_novo
X_diferente <- cbind(const = 1, t = t_novo)
X_diferente

#número de condição da X_padronizada:
kappa(X_diferente) # muito melhor!

# Solução de quadrados mínimos via QR:
#decomposição QR:
QR_d <- qr(X_diferente)
#pegando a Q:
Q_d <- qr.Q(QR_d)
Q_d
#pegando a R:
R_d <- qr.R(QR_d)
R_d

# resolver:   X = Q R 
#             X beta = y
#           (Qtransp X) beta = Qtransp y
#             R beta = Qtransp y
# solucao.betas.d <- backsolve(R_d, t(Q_d) %*% y)  ou 
solucao.betas.d <- qr.coef(QR_d, y) 

solucao.betas.d

#Resíduo e sua norma 2:
res_d <- y - X_diferente %*% solucao.betas.d
res_d
res_norm_d <- norm(res_d, "2")
res_norm_d

#------------------------------------------------------------------------------------------
#e)Obtain the normal equations ÂT Âˆx = Ât ˆb. Find the condition number of
#the coefficient matrix. Solve the normal equations. Calculate the norm of the
#residual (\\b — Ax\\2).
#------------------------------------------------------------------------------------------

X_diferente
XdTXd <- t(X_diferente) %*% X_diferente
XdTXd
XdTY <- t(X_diferente) %*% y
XdTY

# número de condição de XTX
kappa(X_diferente)
kappa(XdTXd)

# solução minimos quadrados via equações normais
# resolver:   Xtransp X beta = Xtransp y 
beta_eqnormais_d <- solve(XdTXd, XdTY)
beta_eqnormais_d

# Norma 2 do resíduo
res_eqnormais_d <- y - X %*% beta_eqnormais_d
res_eqnormais_d_norm <- norm(res_eqnormais_d, type = "2")


#------------------------------------------------------------------------------------------
#f)Compare your solutions from parts (d) and (e). This time they are identical (to
#15 decimal places if you care to check). What polynomial do these solutions
#represent? You might want to plot it just to make sure you've worked the
#problem correctly. If you do make a plot, you will see that your line looks
#identical to the one computed using the QR decomposition in part (a). That one
#was correct to only about eight decimal places, whereas the solutions computed
#using the well-conditioned matrices are correct to machine precision, but of
#course the difference does not show up in a plot (and wouldn't matter in most
#applications either).
#------------------------------------------------------------------------------------------

solucao.betas.d
beta_eqnormais_d

comparacao_d <- data.frame(
  metodo = c("QR", "NormalEQ"),
  beta0 = c(solucao.betas.d[1], beta_eqnormais_d["const", 1]),
  beta1 = c(solucao.betas.d[2], beta_eqnormais_d["t", 1]),
  row.names = NULL)
comparacao_d


beta0_qr_d <- comparacao_d$beta0[1]
beta1_qr_d <- comparacao_d$beta1[1]

beta0_ne_d <- comparacao_d$beta0[2]
beta1_ne_d <- comparacao_d$beta1[2]

library(ggplot2)
library(dplyr)

data

df_lines_d <- data.frame(
  t_novo = t_novo,
  y = y,
  y_pred = c(beta0_qr_d + beta1_qr_d * t_novo, beta0_ne_d + beta1_ne_d * t_novo),
  metodo = rep(c("QR", "Normal Equations"), each = length(t_novo))
)
df_lines_d

grafico2 <- df_lines_d %>% ggplot() +
  geom_point(aes(x = t_novo, y = y), color = "black") +
  geom_line(aes(x = t_novo, y = y_pred, color = metodo, linetype = metodo), size = 1.2) +
  scale_color_manual(values = c("QR" = "blue", "Normal Equations" = "red")) +
  scale_linetype_manual(values = c("QR" = "solid", "Normal Equations" = "dashed")) +
  labs(
    title = "Comparação: ajuste Eq lineares  QR vs Equações Normais",
    x = "t_novo",
    y = "y"
  ) +
  theme_minimal(base_size = 14)
grafico2

library(patchwork)
grafico1+grafico2


cat("Sistema (base 1, t) - QR:\n"); 
cat("beta0 =", format(beta0_qr, digits = 15), "\n")
cat("beta1 =", format(beta1_qr, digits = 15), "\n")
cat("\nNorma do resíduo QR:", format(res_norm, digits = 15), "\n")
cat("Sistema (base 1, t) - NormalEq:\n"); 
cat("beta0 =", format(beta0_ne, digits = 15), "\n")
cat("beta1 =", format(beta1_ne, digits = 15), "\n")
cat("\nNorma do resíduo NormalEq:", format(res_eqnormais_norm, digits = 15), "\n\n")


cat("Sistema convertidos para (base 1, 3*10^7*(t-0.98765435)) - QR:\n")
cat("beta0 =", format(beta0_qr_d, digits = 15), "\n")
cat("beta1 =", format(beta1_qr_d, digits = 15), "\n\n")
cat("\nNorma do resíduo QR:", format(res_norm_d, digits = 15), "\n")
cat("Sistema convertidos para (base 1, 3*10^7*(t-0.98765435)) - NormalEq:\n")
cat("beta0 =", format(beta0_ne_d, digits = 15), "\n")
cat("beta1 =", format(beta1_ne_d, digits = 15), "\n\n")
cat("\nNorma do resíduo NormalEq:", format(res_eqnormais_d_norm, digits = 15), "\n\n")
















