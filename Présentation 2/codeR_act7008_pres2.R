### 
### Olico code
###

### ---------------
## Main function
### --------------

estimation_VaR_TVaR <- function(data_list, kappa, alpha){
  ## calculate VaR for each  sample
  VaR <- lapply(data_list, function(dat) quantile(dat, kappa, type = 4))
  mean_VaR <- mean(unlist(VaR)) # mean VaR
  sd_VaR <- sd(unlist(VaR))
  
  ## Calculate the TVaR for each sample given the VaR
  TVaR <- lapply(1:length(data_list), function(i){
    dat <- data_list[[i]]
    vaR <- VaR[[i]]
    mean(dat[dat > vaR])
  })
  mean_TVaR <- mean(unlist(TVaR))
  sd_TVaR <- sd(unlist(TVaR))
  
  to_return <- data.frame(mean_VaR + sd_VaR / sqrt(length(data_list)) *
                            qt(c(0.025, 0.5, 0.975),
                              df = length(data_list) - 1) ,
             mean_TVaR + sd_TVaR/ sqrt(length(data_list)) *
               qt(c(0.025, 0.5, 0.975),
                  df = length(data_list) - 1))
  
  names(to_return) <- c(paste0("VaR_", kappa*100), paste0("TVaR_", kappa*100))
  rownames(to_return) <- c(paste0(alpha*100, "_lowerbound"),
                           "best estimate",
                           paste0(alpha*100, "_upperbound"))
  return(
         list("table" = to_return,
              "granular_results" = list("VaR" = unlist(VaR),
                                        "TVaR" = unlist(TVaR)))
         )
  
}







### ---------------
## Exemple Oli
### --------------

m <- 1e2   # Nombre de groupe de données
n <- 1e4   # Nombre de données par groupe
kap <-  0.99   # Quantile pour lequel on veut estimer nos VaR/TVaR
alph <- 0.95   # Niveau de confiance désiré pour l'estimation

set.seed(8234)

## simulate the data from a given law
## Je vais faire une binomiale négative composée exponentielle (mélange d'erlang)
data_oli <- lapply(1:m, function(i){
  rgamma(n, rnbinom(n, size = 4, mu = 5), 0.1)
})

## estimation
example_oli <- estimation_VaR_TVaR(data_list = data_oli,
                    kappa = kap,
                    alpha = alph)

p_dist_oli <- function(x){
  dnbinom(0, size = 4, mu = 5) + sum(dnbinom(1:50, size = 4, mu = 5) *
                                       pgamma(x, 1:50, 0.1))
}

VaR_MixErl <- optimize(function(x){log(abs(kap - p_dist_oli(x)))},
                       interval = c(170, 185))$minimum

TVaR_MixErl <- 1/(1 - kap) * sum(sapply(1:50, function(k){
  dnbinom(k, size = 4, mu = 5) * k/0.1 * pgamma(VaR_MixErl, k + 1,
                                                0.1, lower.tail = FALSE)
}))

rbind(example_oli$table,
      "real" = c(VaR_MixErl, TVaR_MixErl))

library(tidyverse)

## boxpplot
example_oli$granular_results %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = L1, y = value)) + geom_boxplot()

library(tidyverse)
library(latex2exp)
## intervalle de confiance
g <- example_oli$granular_results %>% 
  as.data.frame %>% 
  mutate(type = "Simulation") %>% 
  add_row(TVaR = TVaR_MixErl,
          VaR = VaR_MixErl,
          type = "Vrai") %>% 
  ggplot(aes(x = VaR, y = TVaR, color = type)) +
  geom_point(size = 5, alpha = 0.7) + 
  annotate("rect", xmin = -Inf,
           xmax = Inf,
           ymin = example_oli$table$TVaR_99[1],
           ymax = example_oli$table$TVaR_99[3],
           alpha = 0.1, fill = "red")+
  annotate("rect", xmin = example_oli$table$VaR_99[1],
           xmax = example_oli$table$VaR_99[3],
           ymin = -Inf,
           ymax = Inf,
           alpha = 0.1, fill = "red")+
  annotate("rect", xmin = example_oli$table$VaR_99[1],
           xmax = example_oli$table$VaR_99[3],
           ymin = example_oli$table$TVaR_99[1],
           ymax = example_oli$table$TVaR_99[3],
           color = "red",
           alpha = 0)+
  theme_bw() + 
  scale_color_brewer(palette = "Paired") +
  theme(legend.title = element_blank()) +
  labs(x = TeX("$VaR_{\\kappa}(X)$"),
       y = TeX("$TVaR_{\\kappa}(X)$"),
       title = TeX("Estimations de ($VaR_{\\kappa}(X)$, $TVaR_{\\kappa}(X)$) pour chaque échantillon"),
       subtitle = TeX("$X\\sim MixErl$(\\underline{\\zeta}, \\beta = 0.1), $\\zeta_k = P(N = k)$ où $N\\sim BinNeg(r = 4, \\mu = 5)$ et \\kappa = 0.99"),
       caption = TeX("Pour les simulations :  \\alpha = 0.95, $n = 10000$ et m = 100")) + 
  annotate(geom = "text", x = 174, y = 215,
           label = "Intervalle de confiance \n sur la VaR") +
  annotate(geom = "curve", x = 175.5, y = 216.5,
           xend = example_oli$table$VaR_99[1] - 0.25,
           yend = 218,
          curvature = -.3, arrow = arrow(length = unit(2, "mm"))) + 
  annotate(geom = "text", x = 183, y = 202,
           label = "Intervalle de confiance \n sur la TVaR") +
  annotate(geom = "curve", x = 184.5, y = 203.5,
           xend = 185,
           yend = example_oli$table$TVaR_99[1] - 0.25,
           curvature = .3, arrow = arrow(length = unit(2, "mm")))

ggsave("fig/pres3_exnum_oli.png", g)

library(xtable)
xtable(rbind(example_oli$table,
             "real" = c(VaR_MixErl, TVaR_MixErl)))

### ---------------
## Exemple XXX
### --------------

m <- 1e3   # Nombre de groupe de données
n <- 1e4   # Nombre de données par groupe
kap <-  0.99   # Quantile pour lequel on veut estimer nos VaR/TVaR
alph <- 0.95   # Niveau de confiance désiré pour l'estimation

## simulate the data from a given law
# data_XXX <- lapply(1:m, function(i){
#   rXXX(n, ... )
# })

## estimation
example_XXX <- estimation_VaR_TVaR(data_list = data_XXX,
                                   kappa = kap,
                                   alpha = alph)

example_XXX$table # votre table

library(xtable) # installer au besoin
xtable(example_XXX$table) # copier ceci dans vos slides

### -----------
## Graph VaR
### -----------

## Pour un mélange d'Erlang
zeta <- dpois(0:50, lambda = 6)
beta <- 0.1
kappa <- 0.95

stoploss_mixerl <- function(d){
  k <- 1:50
  sum(zeta[k + 1] * (k/beta * pgamma(d, k + 1, beta,
                                      lower.tail = FALSE) -
                       d * pgamma(d, k, beta,
                                  lower.tail = FALSE)))
}

x <- seq(0, 250, by = 0.1)
term1 <- sapply(x, stoploss_mixerl) / (1 - kappa)
term2 <- x
total <- (term1 + term2)

library(tidyverse)
data.frame(x,
           term1,
           term2,
           total) %>% 
  reshape2::melt(id.vars = "x") %>% 
  ggplot(aes(x = x, y = value, color = variable)) +
  geom_line(size = 2, alpha = 0.6) + 
  scale_color_brewer(palette = "Dark2") + 
  # ylim(0, 30)+
  # xlim(40, 175) +
  theme_bw()  
  annotate(geom = "curve", x = 175.5, y = 216.5,
         xend = example_oli$table$VaR_99[1] - 0.25,
         yend = 218,
         curvature = -.3, arrow = arrow(length = unit(2, "mm"))) + 
  annotate(geom = "text", x = 183, y = 202,
           label = "Intervalle de confiance \n sur la TVaR") +
