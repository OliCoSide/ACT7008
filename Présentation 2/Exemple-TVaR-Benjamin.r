##Exemple : somme de deux Pareto
#
library(actuar)

m <- 1000000

set.seed(20)

U1 <- runif(m)
U2 <- runif(m)

X1 <- qpareto(U1, shape = 2, scale = 100)
X2 <- qpareto(U2, shape = 3, scale = 200)

SS <- X1+X2
SSo <- sort(SS)


##kappa = 0.9
k <- 0.9

##VaR
VaR <- SSo[k*m]
VaR

#intervalle de confiance à 95%
k0 <- floor(qnorm((1+0.95)/2)*sqrt(m*k*(1-k))+0.5)
IC.lower <- SSo[k*m-k0]
IC.upper <- SSo[k*m+k0]

##TVaR
TVaR <- mean(SSo[900001:1000000])
TVaR

#intervalle de confiance à 95%
Variance <- (mean(pmax(SSo-VaR,0)^2)-mean(pmax(SSo-VaR,0))^2)/((1-k)^2*m)
IC.lower <- TVaR-qnorm((1+0.95)/2)*sqrt(Variance)
IC.upper <- TVaR+qnorm((1+0.95)/2)*sqrt(Variance)


##kappa = 0.99
k <- 0.99

##VaR
VaR <- SSo[k*m]
VaR

#intervalle de confiance à 95%
k0 <- floor(qnorm((1+0.95)/2)*sqrt(m*k*(1-k))+0.5)
IC.lower <- SSo[k*m-k0]
IC.upper <- SSo[k*m+k0]

##TVaR
TVaR <- mean(SSo[990001:1000000])
TVaR

#intervalle de confiance à 95%
Variance <- (mean(pmax(SSo-VaR,0)^2)-mean(pmax(SSo-VaR,0))^2)/((1-k)^2*m)
IC.lower <- TVaR-qnorm((1+0.95)/2)*sqrt(Variance)
IC.upper <- TVaR+qnorm((1+0.95)/2)*sqrt(Variance)








