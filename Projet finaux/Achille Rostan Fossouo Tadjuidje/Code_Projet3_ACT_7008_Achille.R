###############################################
## CODE PROJET 3 ACT 7008
## SESSION : HIVER 2022
##
## Auteur : Achille Rostan Fossouo Tadjuidje
###############################################

library(actuar)  
library(tidyverse)
library(stringr)
library(latex2exp)
library(xtable)

###
### Question 1 : Distribution Poisson Composée multivariée
###

to_n <- c(1, 10, 100, 1000)
to_gamma <- c(0.3, 0.2, 0.1, 0.15, 0.25)
to_kappa <- c(0.01, 0.5, 0.99)
lambda <- 1/tan(pi / 8)
to_alpha <- c(0, 0.5, 0.99)
beta <- 0.1
to_ro <- c(0.01, 0.05)

## Espérance de Wn 
 
coefs_v <- function(n, alpha, lam = lambda, to_gam = to_gamma, m = 2^8) {
  
  P.C <- function(s) {
    a1 <- (n * (1 - alpha) + 2 * alpha) / (n + alpha)
    a2 <- alpha * (n - 1) / (n + alpha)
    
    a1 * s + a2 * s^2
  }
  #P.C(0.5, 10, 0.3)
  
  P.Nn <- function(s) {
    lam1 <- lam * (n + alpha) / (1 + alpha) 
    exp(lam1 * (P.C(s) - 1))
  }
  
  #P.Nn(0.5, 10, 0.3)
  
  P.J <- function(s) {
    to_sum <- sapply(1:length(to_gam), function(j) {
      to_gam[j] * s^j
    })
    sum(to_sum)
  }
  
  P.K <- function(s) { 
    P.Nn(P.J(s))
  }
  
  # ft_Wn <- function(t) {
  #   s <- n * bet / (n * bet + t)
  #   P.K(s)
  # }
  
  P.K <- Vectorize(P.K)
  
  deg <- c(0, 1) # sévérité dégénére pour avoir la fréquence
  deg_long <- c(deg, rep(0, m - length(deg))) # On ajoute des 0
  fft_deg <- fft(deg_long) # fonction caractéristique de deg
  
  fft_v <- P.K(s = fft_deg)
  
  Re(fft(fft_v, inverse = TRUE))/m   
}

E.Wn <- function(n, alpha, bet = beta) {
  bet2 <- bet * n
  vj <- coefs_v(n, alpha)
  sum(vj * (0:(length(vj) - 1))) / bet2
}

table_E.Wn <- sapply(to_n, function(tn) {
  sapply(to_alpha, function(al) {
    E.Wn(tn, al)
  })
})

colnames(table_E.Wn) <- paste("n = ", to_n)
rownames(table_E.Wn) <- paste("alpha = ", to_alpha)

Var.Wn <- function(n, alpha, bet = beta) {
 #E.Wn.k(n, alpha, k = 2) - E.Wn.k(n, alpha, k = 1)^2
  bet2 <- bet * n
  vj <- coefs_v(n, alpha)
  ek <- sum((0:(length(vj) - 1)) * vj)
  vk <- sum(((0:(length(vj) - 1)) - ek)^2 * vj)
  
  (vk + ek) / bet2^2
}

table_Var.Wn <- sapply(to_n, function(tn) {
  sapply(to_alpha, function(al) {
    Var.Wn(tn, al)
  })
}) 

colnames(table_Var.Wn) <- paste("n = ", to_n)
rownames(table_Var.Wn) <- paste("alpha = ", to_alpha)

## Plus simplement...
## utiliser les questions 7 et 8

## Espérance de Wn
# E(Wn) = lambda * E(B)
E.B <- sum(to_gamma * (1:5)/ beta)
(E.Wn <- E.B * lambda) 

## Variance de Wn
E.J <- sum(to_gamma * (1:5))
Var.J <- sum((1:5 - E.J)^2 * to_gamma)

Var.B <- (Var.J + E.J) / beta^2
E.B2 <- Var.B + E.B^2

g.n <- function(n, alpha, lam = lambda) {
  (lam * E.B2 / n) + (2 * (n - 1) / n^2) * (alpha * lam / (1 + alpha)) * E.B^2 
}

table2_Var.Wn <- sapply(to_n, function(tn) {
  sapply(to_alpha, function(al) {
    g.n(tn, al)
  })
})

colnames(table2_Var.Wn) <- paste("n = ", to_n)
rownames(table2_Var.Wn) <- paste("alpha = ", to_alpha)

#xtable(table2_Var.Wn, digits = 4)

## Mesure entropique de Wn
psi.Wn <- function(rho, n, bet = beta, to_gam = to_gamma) {
  s <- n * bet / (n * bet - rho)
  to_sum <- sum(to_gam * (s^(1:5)))
  
  log(to_sum) / rho
}

table_psi.Wn <- sapply(to_n, function(tn) {
  sapply(to_ro, function(rho) {
    psi.Wn(rho, tn)
  })
})

colnames(table_psi.Wn) <- paste("n = ", to_n)
rownames(table_psi.Wn) <- paste("rho = ", to_ro)

#xtable(table_psi.Wn, digits = 4)

## Répartition de MixErl
p_dist <- function(x, n, alpha, bet = beta) {
  bet2 <- bet * n
  vj <- coefs_v(n, alpha) 
  len <- length(vj)
  vj[1] + sum(vj[2:len] * pgamma(x, 1:(len - 1), bet2))
}

VaR_MixErl <- function(kap, n, alpha, interval = c(0, 69)) {
  optimize(function(x) {
    log(abs(kap - p_dist(x, n, alpha)))
  }, interval = interval)$minimum
}

VaR_MixErl(0.01, 1, 0)

table_VaR_MixErl <- function(kap) {
  
  t1 <- sapply(to_n, function(tn) {
    sapply(to_alpha, function(al) {
      VaR_MixErl(kap, tn, al)
    })
  }) 
  
  colnames(t1) <- paste("n = ", to_n)
  rownames(t1) <- paste("alpha = ", to_alpha)
  t1
}

tab1 <- table_VaR_MixErl(0.01)
tab2 <- table_VaR_MixErl(0.5)
tab3 <- table_VaR_MixErl(0.99)


#xtable(tab1, digits = 4)
#xtable(tab2, digits = 4)

TVaR_MixErl <- function(kap, n, alpha, bet = beta) {
  bet2 <- n * beta
  vj <- coefs_v(n, alpha) 
  len <- length(vj)
  
  varf <- VaR_MixErl(kap, n, alpha)
  
  ss1 <-  sum(vj[2:len] * (1:(len - 1)) * (1 - pgamma(varf, 2:len, bet2)) / bet2)
 
  sum(ss1) / (1 - kap)
}

TVaR_MixErl(0.01, 10, 0.5)

table_TVaR_MixErl <- function(kap) {
  t1 <- sapply(seq_along(to_n), function(tn) {
    sapply(seq_along(to_alpha), function(al) {
      TVaR_MixErl(kap, tn, al)
    })
  }) 
  
  colnames(t1) <- paste("n = ", to_n)
  rownames(t1) <- paste("alpha = ", to_alpha)
  t1
}

tav1 <- table_TVaR_MixErl(0.01)
tav2 <- table_TVaR_MixErl(0.5)
table_TVaR_MixErl(0.99)

#xtable(tav1, digits = 4)
#xtable(tav2, digits = 4)

###
### Question 2 : Modèle basé sur une copule archimédienne imbriquée
###

alpha_j <- 0.1 * (0:2) + 0.2 #j = 0, 1, 2

#################################
#################################

k_max <- 11

pi_s_Xi <- function(s, i) {
  if (s > 0 & s < 1){
  return(sum((1 - pbinom(0:10, size = 10, prob = 0.05 * i))^s))
  }
  if (s >= 1){
    return(sum(1 - (pbinom(0:10, size = 10, prob = 0.05 * i))^s))
  }
}

vect_s <- c(0.5, 0.8, 1, 1.25, 2)
vect_i <- 1:5

pis <- sapply(vect_s, function(s) {
  sapply(vect_i, function(i) {
    pi_s_Xi(s, i)
  })
})

colnames(pis) <- paste('s = ', vect_s)
rownames(pis) <- paste('i = ', vect_i)
pis

#xtable(pis, digits = 4)

#TLS inverse de Theta
alpha_j <- 0.1 * (0:2) + 0.2 #j = 0, 1, 2

#Fonction de répartition de Xij
F.Xij <- function(i, j, kmax = k_max) {
  # i = 1, j = 1, 2
  # i = 2, j = 1, 2, 3
  
  l <- 1 * (i == 1 & j == 1) + 2 * (i == 1 & j == 2) + 3 * (i == 2 & j == 1) 
  + 4 * (i == 2 & j == 2) + 5 * (i == 2 & j == 3)
  
  pbinom(0:kmax, size = 10, prob = 0.05 * l)
}

TLS.inv <- function(i, u1, alph = alpha_j){
  # i = 1, 2
  al <- alpha_j[i + 1]
  log((1 - al)/u1 + al)
}

e_X1u_X2v <- function(u, v, nmax = 100 ) {
  #u = 1, 2
  #v = 1, 2, 3
  
  f_X1u_th0 <- function(k, u, th0) {
    # u = 1, 2
    tot_to <- function(x) exp(-th0 * TLS.inv(i = 1, u1 = pbinom(x, size = 10, prob = 0.05 * u)))
    if (k == 0) {
      return(tot_to(0))
    }
    if (k > 0) {
      return(tot_to(k) - tot_to(k - 1))
    }
  }
  
  f_X2v_th0 <- function(k, v, th0) {
    # v = 1, 2, 3
    l <- 3 * (v == 1) + 4 * (v == 2) + 5 * (v == 3)
    tot <- function(x) exp(-th0 * TLS.inv(i = 2, u1 = pbinom(x, size = 10, prob = 0.05 * l)))
    if (k == 0) {
      return(tot(0))
    }
    if (k > 0) {
      return(tot(k) - tot(k - 1))
    }
  }
  
  to_th0 <- sapply(1:nmax, function(t0) {
    
    sum(sapply(0:nmax, function(k1u) {
      k1u * f_X1u_th0(k1u, u, t0)
    })) *
    sum(sapply(0:nmax, function(k2v) {
        k2v * f_X2v_th0(k2v, v, t0)
    })) *
    0.2^(t0 - 1) * 0.8  
    
  })
  
  sum(to_th0)
}
e_X1u_X2v(2, 2) 

cov_X1u_X2v <- function(u, v) {
  # u = 1, 2
  # v = 1, 2, 3
 
  ex1u <- 10*0.05*u
  ex2v <- 10*0.05*(v + 2)
  
  e12 <- e_X1u_X2v(u, v)
  
  e12 - ex1u * ex2v
}

e_X11_X12 <- function(nmax = 100) {
 
  f_X11_th1 <- function(k, th1) { 
    tot_to <- function(x) exp(-th1 * TLS.inv(i = 1, u1 = pbinom(x, size = 10, prob = 0.05)))
    if (k == 0) {
      return(tot_to(0))
    }
    if (k > 0) {
      return(tot_to(k) - tot_to(k - 1))
    }
  }
  
  f_X12_th1 <- function(k, th1) {  
    tot <- function(x) exp(-th1 * TLS.inv(i = 1, u1 = pbinom(x, size = 10, prob = 0.05 * 2)))
    if (k == 0) {
      return(tot(0))
    }
    if (k > 0) {
      return(tot(k) - tot(k - 1))
    }
  }
  
  to_th1 <- sapply(1:nmax, function(t1) {
    
    sum(sapply(0:nmax, function(k11) {
      k11 * f_X11_th1(k11, t1)
    })) *
      sum(sapply(0:nmax, function(k12) {
        k12 * f_X12_th1(k12, t1)
      })) *
      0.3^(t1 - 1) * 0.7  
    
  })
  
  sum(to_th1)
}

cov_X11_X12 <- e_X11_X12() - 0.5

e_X2u_X2v <- function(u, v, nmax = 100 ) {
  #u = 1, 2, 3
  #v = 1, 2, 3
  # u différent de v
  
  f_X2u_th2 <- function(k, u, th2) {
    # u = 1, 2, 3
    tot_to <- function(x) exp(-th2 * TLS.inv(i = 2, u1 = pbinom(x, size = 10, prob = 0.05 * (u + 2))))
    if (k == 0) {
      return(tot_to(0))
    }
    if (k > 0) {
      return(tot_to(k) - tot_to(k - 1))
    }
  }
  
  f_X2v_th2 <- function(k, v, th2) {
    # v = 1, 2, 3 
    tot <- function(x) exp(-th2 * TLS.inv(i = 2, u1 = pbinom(x, size = 10, prob = 0.05  * (v + 2))))
    if (k == 0) {
      return(tot(0))
    }
    if (k > 0) {
      return(tot(k) - tot(k - 1))
    }
  }
  
  to_th2 <- sapply(1:nmax, function(t2) {
    
    sum(sapply(0:nmax, function(k2u) {
      k2u * f_X2u_th2(k2u, u, t2)
    })) *
      sum(sapply(0:nmax, function(k2v) {
        k2v * f_X2v_th2(k2v, v, t2)
      })) *
      0.4^(t2 - 1) * 0.6  
    
  })
  
  sum(to_th2)
}

e_X2u_X2v(1, 2)

cov_X2u_X2v <- function(u, v) {
  e22 <- e_X2u_X2v(u, v)
  # u = 1, 2, 3
  # v = 1, 2, 3
  
  ex2u <- 10*0.05*(v + 2)
  ex2v <- 10*0.05*(v + 2) 
  
  e22 - ex2u * ex2v
}

matc <- matrix(NA, nrow = 5, ncol = 5)

for (i in 1:5) {
  for (j in i:5) {
    if (i == j) {
      matc[i, j] <- 0.5 * i * (1 - 0.05 * i)
    }
  }
}

matc[1, 2] <- cov_X11_X12
matc[1, 3:5] <- cov_X1u_X2v(1, 1:3)
matc[2, 3:5] <- cov_X1u_X2v(2, 1:3)
matc[3, 4:5] <- cov_X2u_X2v(1, 2:3)
matc[4, 5] <- cov_X2u_X2v(2, 3)

matc[lower.tri(matc)] <- t(matc)[lower.tri(matc)]
xtable(matc, digits = 4)

sum(matc) - sum(sapply(1:5, function(i) matc[i, i]))

cov_xi_s <- function(i) {
  vi <- matc[i, i]
  sci <- 0
  for (j in 1:5) {
    if (j != i) sci <- sci + matc[i, j]
  }
  vi + sci
}

vs <- 41.42771

ft <- sapply(1:5, function(i) cov_xi_s(i))
sum(ft / sqrt(vs))
