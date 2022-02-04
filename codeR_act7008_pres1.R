##
## Olivier
## 

## packages ---------

library(tidyverse)


## fonctions --------------
## Dekingkong
c_dk <- function(j, k, n, vs, vf, a, b){
  term1 <- vs * choose(n - 2 - j, k - 1)
  term2 <- (vs * a + vf * b)/(1 - b) * choose(n - 2 - j, k)
  term3 <- (vf * a * b)/(1 - b)^2 * choose(n - 2 - j, k + 1)
  return(term1 + term2 + term3)
}

## Dekingkong
del <-function(a, b){
  return(
    a*b/((1 - a) * (1 - b))
  )
}

## Dekingkong
f_M <- function(n, j, vs, vf = 1 - vs, a, b){
  ## Off target
  if(j <0 || j > n)return(0)
  
  ## min value
  if(j == 0){
    return(vf * (1 -b)^(n - 1))
  }
  
  ## max value
  if(j == n){
    return(vs * (1 -a)^(n - 1))
  }
  
  terms <- sapply(0:(j - 1), function(k){
    choose(j - 1, k) * del(a, b)^k * c_dk(j = j-1, k = k,
                                       n = n, vs = vs, vf = vf,
                                       a= a, b = b)
  })
  
  return(
    (1 - b)^(n - j) * (1 - a)^(j - 1) * sum(terms)
  )
}

## Notation Cossette 2003 dans la fmp de Dekingkong
f_M_applied <- function(n, j_tot, p, pi){
  
  mp <- matrix_bm(p = p, pi = pi, s = 1,
                            power = 1, eig = FALSE)
  
  vec_dens <- sapply(j_tot, function(j){
    term1 <- ifelse(j == 0, mp[1,1]^n,
                    ifelse(j == n,
                           mp[1, 2] * mp[2,2]^(n - 1),
                           f_M(n, n - j, vs = mp[1,1],
                               vf = mp[1, 2],
                               a = (1 - pi) * p,
                               b = (1 - pi) * (1 - p))))
    term2 <- ifelse(j == 0, mp[2, 1] * mp[1,1]^(n - 1),
                    ifelse(j == n,
                           mp[2,2]^n,
                           f_M(n, n - j, vs = mp[2,1],
                               vf = mp[2, 2],
                               a = (1 - pi) * p,
                               b = (1 - pi) * (1 - p))))
    (1 - p) * term1 + p * term2
  })
  
  return(
    vec_dens
  )
}

## fgp Cossette et collab 2003 avec la formule de Deking kong
P_M_expl <- function(n, s_tot, p, pi){
  sapply(s_tot, function(s){
    probs <- sapply(0:n, function(j) f_M_applied(n, j, p = p, pi = pi))
    sum(s^(0:n) * probs)
  })
}

matrix_bm <- function(p, pi, s, power = 1,
                      eig = FALSE){
  
  if(power <= 0) return(diag(2))
  
  p_00 <- 1- (1-pi) * p
  p_01 <- (1 - pi) * p * s
  p_10 <- (1 - pi) * (1 - p)
  p_11 <- ((1 - pi) * p + pi)  * s
  
  ## initial matrix
  mat <- matrix(c(p_00, p_01,
                  p_10, p_11),
                byrow = TRUE,
                nrow = 2)
  
  ## if eigen
  if(eig == TRUE){
    a <- eigen(mat)
    return(a$vectors %*% diag(a$values^(power)) %*% solve(a$vectors))
  }
  
  mat_to_return <- mat %*% matrix_bm(p, pi, s,
                                     power = power - 1,
                                     eig = FALSE)
  return(mat_to_return)
}

## fgp de
P_M_mat <- function(n, s_tot, p, pi, eig = FALSE){
  sapply(s_tot, function(s){
    left <- matrix(c(1 - p, p * s), nrow = 1)
    middle <- matrix_bm(p, pi, s, power = n - 1,
                        eig = eig)
    right <- matrix(c(1, 1), nrow = 2)
    (left %*% middle %*% right)[1, 1]
  })
}

P_Z <- function(s, pi){
  (1 - pi) * s / (1 - pi * s)
}

P_K <- function(s, lam){
  exp(lam * (s - 1))
}

P_N <- function(n, s_tot, p, pi){
  lam <- n * p
  sapply(s_tot, function(ss){
    P_K(P_Z(ss, pi), lam = (1 - pi) * lam)
  })
}

## test ------

pi_to_graph <- c(0, 0.2, 0.4, 0.6, 0.8)
names(pi_to_graph) <- sapply(pi_to_graph *100, function(p) paste0("ex_", p))

table1_mk <- sapply(pi_to_graph, function(piii){
  ## tableau 1 cossette 2003
  f_const <- c(0, 1)
  aa <- 2^8
  nb <- length(f_const) # length right now
  ftc <- fft(c(f_const, rep(0, aa - nb))) # fonction caractéristique
  f_S_n <- Re(fft(P_M_mat(20, ftc, 0.1, pi = piii, eig = TRUE), TRUE))/aa # On inver
  return(round(f_S_n, 6)[1:21])
})

## on veut la table LaTeX
library(xtable)
xtable(table1_mk[9:12, 1:5], digits = 6)


table1_mk2 <- sapply(pi_to_graph, function(piii){
  ## tableau 1 cossette 2003
  f_const <- c(0, 1)
  aa <- 2^8
  nb <- length(f_const) # length right now
  ftc <- fft(c(f_const, rep(0, aa - nb))) # fonction caractéristique
  f_S_n <- Re(fft(P_M_expl(20, ftc, 0.1, pi = piii), TRUE))/aa # On inver
  return(round(f_S_n, 6)[1:21])
})

## MIeux
table1_mk2 <- sapply(pi_to_graph, function(piii){
  ## Nous avons le format explicite
  f_M_k <- f_M_applied(20, 1:1e3, p = 0.1,
                       pi = piii) 
  return(round(f_M_k, 6)[1:21])
})

pi_to_graph2 <- c(0, 0.2, 0.4, 0.6, 0.8)
names(pi_to_graph2) <- sapply(pi_to_graph *100, function(p) paste0("ap_", p))


table1_mk3 <- sapply(pi_to_graph2, function(piii){
  ## tableau 1 cossette 2003
  f_const <- c(0, 1)
  aa <- 2^8
  nb <- length(f_const) # length right now
  ftc <- fft(c(f_const, rep(0, aa - nb))) # fonction caractéristique
  f_S_n <- Re(fft(P_N(20, ftc, 0.1, pi = piii), TRUE))/aa # On inver
  return(round(f_S_n, 6)[1:21])
})

values_x_comp <- seq(from = 0, by = 1, length.out = nrow(table1_mk3))
## comparison, pois-non pois
data_comparison <- data.frame("x" = values_x_comp,
                              table1_mk3,
                              table1_mk)

data_comparison2 <- data_comparison %>% reshape2::melt(id.vars = "x") %>% 
  mutate(pi = (variable %>% substr(4, 5) %>% as.numeric)/100,
         type = variable %>% substr(1, 2)) %>% 
  select(-variable)

data_comparison2 %>% ggplot(aes(x = x, y = value,
                                     fill = factor(pi))) +
  geom_bar(stat= "identity", position = "dodge")

ggsave("plot/approx.png",
       data_comparison2 %>%  
         filter(pi == 0.4) %>% 
         ggplot(aes(x = x,
                    color = type,
                    fill = type,
                    y = value)) + 
         geom_bar(stat = "identity", 
                  position = "dodge",
                  lwd = 0, alpha = 0.6) +
         theme_bw() + 
         scale_color_brewer(palette = "Dark2", labels = c("Approximation", "Exact")) +
         
         scale_fill_brewer(palette = "Dark2", labels = c("Approximation", "Exact")) +
         labs(x = "s",
              y = "Densité",
              title = TeX("Comparaison entre $f_{M_{20}}(s)$ et $f_{N}(s)$"),
              subtitle = TeX("Modèle binomial markovien avec $\\pi = 0.4$ et $q = 0.1$")) + 
         scale_y_continuous(labels = scales::percent) + 
         theme(legend.title = element_blank())
)

## Soit B qui suit une binomiale négative (r = 5, q = 0.1)
f_M_applied(20, 0:5, p = 0.1, pi = 0.4) %>%  round(6)


eps <- 1e-10

## 1. trouver une valeur b* tq. F_B >= 1-1e-10
bmax <- qnbinom(1-eps, 5, 0.1) + 1

## densité B
f_B <- sapply(1:bmax, function(x) dnbinom(x - 1, 5, 0.1))


p <- 0.1
pi <- 0.3
n <- 15

aa <- 2^18
nb <- length(f_B) # length right now
ftb <- fft(c(f_B, rep(0, aa - nb))) # fonction caractéristique
f_S_n <- Re(fft(P_M_mat(n, ftb, p, pi, eig = TRUE), TRUE))/aa # On inver

library(actuar)

## Ex cossette 2003
beta <- 26.519019
f_log <- sapply(1:2e2, function(j) beta^j/(j * (1 + beta)^j * log(1 + beta)))

pi_for_table2_basic <- c(0, 0.4, 0.8)
names(pi_for_table2_basic) <- sapply(pi_for_table2_basic *100, function(p) paste0("pi_", p))

table2basic <- sapply(pi_for_table2_basic, function(piiii){
  aa <- 2^10
  nb <- length(f_log) # length right now
  ft_log <- fft(c(0, f_log, rep(0, aa - 1 - nb))) # fonction caractéristique
  f_S_n <- Re(fft(P_M_mat(100, ft_log, p = 0.1, pi = piiii, eig = TRUE), TRUE))/aa # On inver
  F_S_n <- cumsum(f_S_n)
  F_S_n[1:301]
})

xtable(table2basic, digits = 6)

pi_for_table2 <- seq(0, 0.98, by = 0.02)
names(pi_for_table2) <- sapply(pi_for_table2 *100, function(p) paste0("pi_", p))

table2 <- sapply(pi_for_table2, function(piiii){
  aa <- 2^10
  nb <- length(f_log) # length right now
  ft_log <- fft(c(0, f_log, rep(0, aa - 1 - nb))) # fonction caractéristique
  f_S_n <- Re(fft(P_M_mat(100, ft_log, p = 0.1, pi = piiii, eig = TRUE), TRUE))/aa # On inver
  F_S_n <- cumsum(f_S_n)
  F_S_n
})

values_x2 <- seq(from = 0, by = 1, length.out = nrow(table2))

dat <- data.frame("x" = values_x2,
                  table2)

dat2 <- dat %>% reshape2::melt(id.vars = "x") %>% 
  mutate(pi = (variable %>% substr(4, 5) %>% as.numeric)/100) 

setwd("C:/Users/olico/OneDrive - Université Laval/ULaval/10 - H22/ACT-7008 Sujets spéciaux III/pres1")

col_vi <- hcl.colors(5, "Plasma")
library(latex2exp)
ggsave("plot/repart.png",
       dat2 %>%  ggplot(aes(x = x,
                            color = pi,
                            group = factor(variable))) + 
         geom_line(aes(y = value), lwd = 2, alpha = 0.6) +
         scale_colour_gradient(name = TeX("Valeur de $\\pi$"),
                               low = col_vi[1],
                               high = col_vi[4],
                               trans = "exp") + 
         theme_bw() + 
         labs(x = "s",
              y = "Fonction de répartition",
              title = TeX("$F_{S_{100}}(s)$ selon la valeur $\\pi$"),
              subtitle = TeX("$q = 0.1$ et sévérité $X \\sim logarithimque(\\beta = 26.519019)$")) + 
         scale_y_continuous(labels = scales::percent) + 
         scale_x_continuous(labels = scales::dollar)
)

library(tictoc)
tic("mat")
f_S_n <- Re(fft(P_M_mat(n, ftb, p, pi, eig = FALSE), TRUE))/aa # On inver
toc()
tic("eig")
f_S_n <- Re(fft(P_M_mat(n, ftb, p, pi, eig = TRUE), TRUE))/aa # On inver
toc()

## general param
k_per <- 50
n_simul <- 50
eta <- 0.15
u <- 13

## param BM
p <- 0.3
pi <- 0.2
lam <- (p*(1 + eta))

simuls <- sapply(1:n_simul, function(i){
  uu <- c(u, rep(0, k_per))
  ii <- c(rbinom(1, 1, p), rep(0, k_per))
  mat <- matrix_bm(p, pi, s =1)
  
  for (k in 1:k_per) {
    prob <- ifelse(ii[k] == 1, mat[2, 2], mat[1, 2])
    
    ii[k + 1] <- rbinom(1, 1, prob)
    b <- rexp(1, lam)
    y <- ii[k + 1] * b
    uu[k + 1] <- uu[k] + 1 - y
  }
  return(uu)
})

periods <- 0:k_per

mean(simuls %>% apply(2, function(t) quantile(t, 0.01)))

data_to_graph <- data.frame(simuls, "period" = periods) %>% 
  reshape2::melt(id.vars = "period") %>% 
  group_by(variable) %>% 
  mutate(quant = quantile(value, 0.01),
         is_ruined = (min(value) < 0))

ggsave("plot/simuls.png",
data_to_graph %>%
  ggplot(aes(x = period, y = value,
             color = quant, group = variable)) +
  geom_line(alpha = 0.5, lwd = 1) + 
  scale_color_gradient(low = "deeppink3",
                       high = "palegreen3") + 
  
  theme_bw() + 
  labs(x = "Période",
       y = "Surplus",
       title = TeX("Simulations de $U_k$ avec modèle binomial markovien composé"),
       subtitle = TeX("$q = 0.3$, $\\pi = 0.2$, $u = 13$, $\\eta = 0.15$ et sévérité $X \\sim Exp(\\lambda = q (1 + \\eta))$")) + 
  scale_y_continuous(labels = scales::dollar) + 
  theme(legend.position = "")
)
#        dat2 %>%  ggplot(aes(x = x,
#                             color = pi,
#                             group = factor(variable))) + 
#          geom_line(aes(y = value), lwd = 2, alpha = 0.6) +
#          scale_colour_gradient(name = TeX("Valeur de $\\pi$"),
#                                low = "deeppink3",
#                                high = "greenyellow",
#                                trans = "exp") + 
#          theme_bw() + 
#          labs(x = "s",
#               y = "Fonction de répartition",
#               title = TeX("$F_{S_k}(s)$ selon la valeur $\\pi$"),
#               subtitle = "Modèle binomial markovien") + 
#          scale_y_continuous(labels = scales::percent) + 
#          scale_x_continuous(labels = scales::dollar)
# )
