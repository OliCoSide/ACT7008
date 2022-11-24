##===================================================
## PROJET2 ACT 7008 - THEORIE DE LA MESURE DU RISQUE
##
## 
##===================================================

## Monte Carlo method
## Avec la loi exponentielle (beta = 1)
 
VaR_TVaR_kap <- function(kappa, alpha = 0.05, nsim = 2^16, seed = 202203) {
  
  set.seed(seed = seed)
  vect.simul <- rexp(nsim)
  
  ## VaR : estimation
  VaR_estim <- sort(vect.simul)[ceiling(nsim * kappa)]
  
  ## VaR : vraie valeur
  VaR_vrai <- - log(1 - kappa)
  
  ## VaR : intervalle de confiance
  ic_var <- VaR_estim + c(-1, 1) * qnorm(1 - alpha/2) * sqrt(kappa * (1 - kappa) / nsim) / dexp(VaR_vrai)
  
  ## TVaR : vraie valeur
  TVaR_vrai <- VaR_vrai + 1
  
  ## TVaR : estimation
  TVaR_estim <- VaR_estim + 1 / (nsim * (1 - kappa)) * sum(pmax(vect.simul - VaR_estim, 0))
  
  ## TVaR : intervalle de confiance
  ic_tvar <- TVaR_estim + c(-1, 1) * qnorm(1 - alpha/2) * sqrt((2 * exp(- VaR_vrai) - exp(- 2 * VaR_vrai)) / nsim) / (1 - kappa)
  
  to_return <- data.frame(c(ic_var[1], VaR_estim, ic_var[2], VaR_vrai),
                          c(ic_var[2], TVaR_estim, ic_tvar[2], TVaR_vrai))
  
  names(to_return) <- c(paste0("VaR_", kappa*100), paste0("TVaR_", kappa*100))
  rownames(to_return) <- c("Borne inférieure",
                           "Valeur estimée",
                           "Borne supérieure",
                           "Vraie valeur")
   
  return(to_return)
}

VaR_TVaR_kap(0.99)


library(xtable) # installer au besoin
xtable(VaR_TVaR_kap(0.99), digits = 6)
