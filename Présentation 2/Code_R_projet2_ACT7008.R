##===================================================
## PROJET2 ACT 7008 - THEORIE DE LA MESURE DU RISQUE
##
## 
##===================================================

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
  
  return(list("VaR" = list("Vraie valeur" = VaR_vrai,
                           "Valeur estimée" = VaR_estim,
                           "Intervalle de confiance" = ic_var),
              
              "TVaR" = list("Vraie valeur" = TVaR_vrai,
                            "Valeur estimée" = TVaR_estim,
                            "Intervalle de confiance" = ic_tvar)))
}

VaR_TVaR_kap(0.95)
