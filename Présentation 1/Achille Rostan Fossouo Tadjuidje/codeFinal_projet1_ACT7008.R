##=========================================================================
## ACT 7008 
##
## Sujets spéciaux - session d'Hiver 2022 
##
## Lien de l'article :  https://doi.org/10.1080/03461238.2021.1937305 
##
## Reproduire l'exemple numérique de la section 5.1 de l'article 
## avec un processus de Poisson INAR(1)
##=========================================================================


library(actuar)
library(ggplot2)
library(tidyverse)
library(stringr)
library(latex2exp)

## n = nombre de périodes
n <- 4

beta <- 1

## lambda = paramètre de la loi Poisson
lambda <- 1
 
## dependance parameter
vect_alpha <- c(0.25, 0.5, 0.75)
alpha <- vect_alpha[1]

vect_kappa <- c(0.99, 0.995)
kappa <- vect_kappa[1]
 
###===========================================
### Apllication avec Sn
###=========================================== 

## FGP de la Loi de Epsilon : Poisson(lambda)
P_eps <- function(z, lambda) exp(lambda * (z - 1))

## FGP de la loi de  N : Poisson(lambda/(1 - alpha))
P_N <- function(z, lam, alpha) P_eps(z, lambda = lam / (1 - alpha)) 

## Fonction g (Equation numéro 20 dans l'article)
g <- function(n, k, alpha) {
  if (k == n) return(alpha^(n - 1))
  if (!(k %in% 1:n)) return(0)
  
  alpha^(k - 1) * (1 - alpha)
}
 

G <- function(x, n, alpha){
  if(n == 1) return(x)
  terms <- sapply(1:n, function(k)g(n, k, alpha) * x^k)
  sum(terms)
}

P_J <- function(t){
  0.4 * t +0.6 * t^2
}

P_Mstar <- function(t, n, lambda, alpha){
  to_prod <- sapply(2:n, function(i){
    P_eps(
      G(P_J(t), n + 1 - i, alpha),
      lambda)
  })
  P_N(G(P_J(t), n, alpha), lambda, alpha) * prod(to_prod)
}
P_Mstar <- Vectorize(P_Mstar)

## FFT
m <-  6# exposant
aa <- 2^m # longueur de vecteur attendu

f_B_deg <- c(0, 1) # sévérité dégénére pour avoir la fréquence

f_B_deg_long <- c(f_B_deg, rep(0, aa - length(f_B_deg))) # On ajoute des 0

fft_B_deg <- fft(f_B_deg_long) # fonction caractéristique de B deg

# fonction caractéristique de Sn (avec sévirité dégénérée, 
# donc meme chose que fonction caractéristique de Mstar)
fft_M_star <- function(alph) P_Mstar(t = fft_B_deg, alpha = alph, n = n, lambda = lambda)

f_M_star <- function(alph) Re(fft(fft_M_star(alph), inverse = TRUE))/aa # On inverse pour obtenir la fmp de Sn_star (donc la fmp de Mstar puisque B est dégénéré à 1)


E.Sn <- function(alpha) {
  fm <- f_M_star(alpha)
  sum( fm * (0:(length(fm) - 1)))
}  

plot.Esn <- data.frame(seq.alpha = seq(0, 0.999, by = 0.001), seq.E.Sn = sapply(seq(0, 0.999, by = 0.001), function(al) E.Sn(al)))

col_vi <- hcl.colors(5, "Plasma") 

ggsave("graph_esperance.png",
plot.Esn %>%
  ggplot(aes(x = seq.alpha, y = seq.E.Sn, colour = seq.alpha)) +
  geom_line(alpha = 0.5, lwd = 2) + 
  scale_color_gradient(low = col_vi[4], high = col_vi[1], na.value = NA) +  
  theme_bw() +
  #labs(x = "alpha",
  #     y = "Espérance de Sn")  +
  labs(x = TeX("alpha"),
       y = TeX("Espérance de Sn"),
       title = TeX("Espérance de Sn selon la valeur de $\\alpha$"),
       subtitle = TeX("n = 4, $\\lambda = 1$, $\\beta = 1$")) +
  theme(legend.position = "")
)

##===================================== Fonction stop-loss
pi.d <- function(d, alpha) {
  ft <- f_M_star(alpha)
  f2 <- sapply(1:length(ft), function(k) {
    (k / beta) * (1 - pgamma(d, k + 1, beta)) - d * (1 - pgamma(d, k, beta))
  })
  
  sum(ft * f2)
}

seq.alpha <- seq(0, 0.99, by = 0.01)
seq.d <- seq(0, 100, by = 1)

data.pi <- sapply(seq.d, function(d) {
  sapply(seq.alpha, function(al) pi.d(d, al))
})

names(seq.alpha) <- sapply(seq.alpha, function(p) paste0("alpha_", p))
names(seq.d) <- sapply(seq.d, function(p) paste0("d_", p))

colnames(data.pi) <- names(seq.d)
rownames(data.pi) <- names(seq.alpha)

#data.pi[c(1:6, 37:41), 1:5]
data.pi <- t(data.pi)



data.pi2 <- as.data.frame(data.pi)
head(data.pi2)
colnames(data.pi2) <- names(seq.alpha)
data.pi2$d <- seq.d 

data.pi2_long <- reshape2::melt(data.pi2, id.vars = "d")

data.pi2_long$alpha <- as.numeric(str_replace(data.pi2_long$variable, "alpha_", ""))



col4 <- hcl.colors(5, "Batlow")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark",  "Spectral"
ggsave("graph_stoploss1.png",
       data.pi2_long %>% 
         ggplot(aes(x = d,
                    color = alpha,
                    group = factor(alpha))) + 
         geom_line(aes(y = value), lwd = 2, alpha = 0.8) +
         scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                               low = col4[1],
                               high = tail(col4, 1),
                               trans = "exp") + 
         theme_bw() + 
         labs(x = TeX("d"),
              y = TeX("$\\pi_{Sn}(d)$"),
              title = TeX("$\\pi_{Sn}(d)$ selon la valeur de $\\alpha$"),
              subtitle = TeX("n = 4, $\\lambda = 1$, $\\beta = 1$")) + 
         scale_y_continuous(labels = scales::dollar) + 
         scale_x_continuous(labels = scales::dollar)
)


## Fonction de répartition de Sn
F.Sn <- function(x, alpha) {
  pk <- f_M_star(alpha)
  
  sum(pk * pgamma(x, 0:(length(pk) - 1), beta))
}

##=========================== VaR
##===
VaR <- function(alpha, kappa) optimise(function(x) abs(F.Sn(x, alpha) - kappa), c(0, 100))$minimum

seq.alpha <- seq(0, 0.99, by = 0.01)
seq.kappa <- seq(0, 1, by = 0.02)

data.var <- sapply(seq.alpha, function(al) {
  sapply(seq.kappa, function(ka) VaR(al, ka))
})

names(seq.alpha) <- sapply(seq.alpha, function(p) paste0("alpha_", p))
names(seq.kappa) <- sapply(seq.kappa, function(p) paste0("kappa_", p))

colnames(data.var) <- names(seq.alpha)
rownames(data.var) <- names(seq.kappa)

#data.pi[c(1:6, 37:41), 1:5]
# data.var <- t(data.var)

data.var2 <- as.data.frame(data.var)
head(data.var2)
#colnames(data.var2) <- names(seq.alpha)
data.var2$kappa <- seq.kappa 

data.var2_long <- reshape2::melt(data.var2, id.vars = "kappa")

data.var2_long$alpha <- as.numeric(str_replace(data.var2_long$variable, "alpha_", ""))

ggsave("graph_var.png",
       data.var2_long %>% 
         ggplot(aes(x = kappa,
                    color = alpha,
                    group = factor(alpha))) + 
         geom_line(aes(y = value), lwd = 2, alpha = 0.8) +
         scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                               low = col4[1],
                               high = tail(col4, 1),
                               trans = "exp") + 
         theme_bw() + 
         labs(x = TeX("kappa"),
              y = TeX("$VaR_{\\kappa}(S_n)$"),
              title = TeX("$VaR_{\\kappa}(S_n)$ selon la valeur de $\\alpha$"),
              subtitle = TeX("n = 4, $\\lambda = 1$, $\\beta = 1$")) + 
         scale_y_continuous(labels = scales::dollar) + 
         scale_x_continuous(labels = scales::dollar)
)

## Fonction de densité de Sn
f.Sn <- function(x, alpha) {
  pk <- f_M_star(alpha)
  
  sum(pk * dgamma(x, 0:(length(pk) - 1), beta))
}
 
##======================= Prime exponentielle de Sn

# FGM :

fgm.b <- function(t, n, beta) (beta / (beta - t))^n

fgm.Sn <- function(z, alpha) {
  pk <- f_M_star(alpha)
  
  sum(pk[-1] * fgm.b(z, 1:(length(pk) - 1), beta)) 
} 

prime.exp <- function(r, alpha) log(fgm.Sn(r, alpha))/r

prime.exp(0.75, alpha)

seq.alpha <- seq(0, 0.99, by = 0.01)
seq.r <- seq(0.01, 0.6, by = 0.01)

data.exp <- sapply(seq.r, function(r) {
  sapply(seq.alpha, function(al) prime.exp(r, al))
})

names(seq.alpha) <- sapply(seq.alpha, function(p) paste0("alpha_", p))
names(seq.r) <- sapply(seq.r, function(p) paste0("r_", p))

colnames(data.exp) <- names(seq.r)
rownames(data.exp) <- names(seq.alpha)

#data.pi[c(1:6, 37:41), 1:5]
# data.var <- t(data.var)
data.exp2 <- t(data.exp)
data.exp2 <- as.data.frame(data.exp2)
#head(data.var2)
#colnames(data.var2) <- names(seq.alpha)
data.exp2$r <- seq.r 

data.exp2_long <- reshape2::melt(data.exp2, id.vars = "r")

data.exp2_long$alpha <- as.numeric(str_replace(data.exp2_long$variable, "alpha_", ""))

ggsave("graph_prime_exp.png",
       data.exp2_long %>% 
         ggplot(aes(x = r,
                    color = alpha,
                    group = factor(alpha))) + 
         geom_line(aes(y = value), lwd = 2, alpha = 0.8) +
         scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                               low = col4[1],
                               high = tail(col4, 1),
                               trans = "exp") + 
         theme_bw() + 
         labs(x = TeX("Aversion au risque r"),
              y = TeX("Prime exponentielle"),
              title = TeX("prime exponentielle selon la valeur de $\\alpha$"),
              subtitle = TeX("n = 4, $\\lambda = 1$, $\\beta = 1$")) + 
         scale_y_continuous(labels = scales::dollar) + 
         scale_x_continuous(labels = scales::dollar)
)
