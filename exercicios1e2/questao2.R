#-----------------------------------------------------------
# matriz A, 2x2
#-----------------------------------------------------------
A <- matrix(c(8, 1,-2, 1), nrow = 2, byrow = TRUE)
A

#autovalores verdadeiros. Veja que o autovalor dominante é 7.70...
eigen(A)$values
#autovetores verdadeiros
eigen(A)$vectors

#-----------------------------------------------------------
#a)Use inverse iteration with p = 8 and qo = [ 1 1 ] to calculate an eigenvalue
#and eigenvector of A. (On this small problem it may be easiest simply to
#calculate B = (A — 8I)^-1 and perform direct iterations with B.)
#-----------------------------------------------------------
p <- 8. # é o deslocamento próximo ao autovalor dominante
q_inicial <- c(1, 1)
q_inicial <- q_inicial/sqrt(sum(q_inicial^2))  #normalizando o vetor inicial
q_inicial

# Calculando B = A — 8I
B <- A - p * diag(2)
B

# Calculando a inversa de B, ou seja, B_inv
B_inv <- solve(B)
B_inv

# Definindo o número de iterações:
n_iter <- 12

# Matriz para guardar os autovetores
Q1 <- matrix(0, nrow = 2, ncol = n_iter+1)
Q1[,1] <- q_inicial

#Vetor para guardar os autovalores
autovalor_estimado <- numeric(n_iter)

for (k in 1:n_iter) {
  q_new <- B_inv %*% q_inicial
  s_new <- q_new[which.max(abs(q_new))]
  q_inicial <- q_new / s_new
  
  Q1[,k+1] <- q_inicial
  
  autovalor_estimado[k] <- as.numeric(t(q_inicial) %*% A %*% q_inicial) #usando quociente de Rayleigh - teorema 5.3.25
}

#Resultados:
autovalor_estimado
Q1

#Sabemos:
#autovalores verdadeiros
eigen(A)$values
#autovetores verdadeiros
eigen(A)$vectors


#-----------------------------------------------------------
#b)Now that you have calculated an eigenvector v, calculate the ratios
#qj+1 - v / qj-v for j — 0, 1, 2, . . . , to find the observed rate of convergence.
#Solve for the eigenvalues using the characteristic equation or MATLAB and calculate the
#theoretical convergence rate | (lambda1 — 8)/(lambda2 — 8) |. How well does theory agree
#with practice?
#-----------------------------------------------------------
#Pegando o autovetor verdadeiro associado ao autovalor verdadeiro mais próximo de p
v_verdadeiro <- eigen(A)$vectors[, which.min(abs(eigen(A)$values - p))] 
v_verdadeiro

#pegando o s para escalonar
s <- max(abs(v_verdadeiro))
s
v_verdadeiro <- v_verdadeiro / s
v_verdadeiro

# sinal:
if (sign(Q1[1,2]) != sign(v_verdadeiro[1])){
  v_verdadeiro <- - v_verdadeiro
} 

#criando um vetor para erros
vetor_erros <- numeric(n_iter+1)

# calculando norma ||q - v|| e salvando no vetor de erros
for (j in 1:(n_iter+1)) {
  vetor_erros[j] <- sqrt(sum((Q1[,j] - v_verdadeiro)^2))
}

vetor_erros

#calculando as razoes:
razoes_observadas <- vetor_erros[-1] / vetor_erros[-length(vetor_erros)]   # e_{j+1}/e_j
razoes_observadas

# taxa teórica:
lambda_1 <- eigen(A)$values[v_verdadeiro]
lambda_1

lambda_2 <- eigen(A)$values[-v_verdadeiro]
lambda_2

taxa_convergencia <- abs((lambda_1 - p) / (lambda_2 - p))
taxa_convergencia 


#Na teoria, temos uma taxa de convergencia muito pequena. Olhando as razões observadas, elas convergem para o valor teórico. 
#Mas, conforme vai aumentando a precisão (os últimos dígitos), essas razões começaram a ficar instáveis.















