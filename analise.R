library(ggplot2)
library(stringi)
library(stringr)
library(data.table)
setwd("/home/ze/Área de Trabalho/")

dados <- fread("trabalho_industria.csv", dec = ",")

dados <- as.data.table(apply(dados, 2, str_replace_all, pattern = " ", replacement = ""))

dados[,val_bolso1 := as.numeric(bolso1)/(as.numeric(bolso1)+as.numeric(haddad1)+as.numeric(ciro)+as.numeric(outros))]
dados[,val_bolso2 := as.numeric(bolso2)/(as.numeric(bolso2)+as.numeric(haddad2))]


cor(as.numeric(dados$variacao_trabalho_industria), dados$val_bolso1, method = "spearman")
cor(as.numeric(dados$evangelicos), dados$val_bolso1, method = "spearman")
cor(as.numeric(dados$catolicos), dados$val_bolso1, method = "spearman")
cor(as.numeric(dados$bovino_18_08), dados$val_bolso1, method = "spearman")
cor(as.numeric(dados$suino_18_08), dados$val_bolso1, method = "spearman")
cor(as.numeric(dados$galo_18_08), dados$val_bolso1, method = "spearman")
cor(as.numeric(dados$var_soja_areacolhida_18_08), dados$val_bolso1, method = "spearman")

u <- rank(as.numeric(dados$evangelicos))/28
v <- rank(as.numeric(dados$val_bolso1))/28

dados %>% ggplot(aes(x = u, y = v)) + 
  geom_point()


u <- rank(as.numeric(dados$var_soja_areacolhida_18_08))/28
v <- rank(as.numeric(dados$val_bolso1))/28

dados %>% ggplot(aes(x = u, y = v)) + geom_point() 

u <- rank(as.numeric(dados$bovino_18_08))/28
v <- rank(as.numeric(dados$val_bolso1))/28

dados %>% ggplot(aes(x = u, y = v)) + geom_point() 


expose_stan_functions("./rnsp_frank_gauss.stan")

## check if copula density integrates to 1
frank_copula_pdf <- function(uv, theta) {
  exp(log_frank_copula(uv[1], uv[2], theta) )
}
adaptIntegrate(frank_copula_pdf
                  ,lowerLimit = c(0,0)
                  ,upperLimit = c(1,1)
                  ,theta = -25
                  ,doChecking = TRUE
                   )
## $integral 1, $error 9.981324e-06
## $functionEvaluations 15113, $returnCode 0

normal.fit.gauss.prior <- stan(
  file = "/home/ze/análise Bolsonaro 2018/normalcopula.stan",
  data = list(
  N = nrow(dados),
  u = qnorm(u),
  v = qnorm(v)
),
seed = 2020
)




