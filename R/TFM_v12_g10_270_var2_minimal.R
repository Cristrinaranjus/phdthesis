

############################################ TÍTULO ############################################
#
# SCRIPT:  SIMULACIÓN TFM - "Aproximaciones bayesianas versus frecuentistas en diseños de caso único desde la perspectiva multinivel:sobre la potencia estadística y la tasa de error tipo I"
# AUTORES: Cristina Rodríguez-Prada y Ricardo Olmos
# FECHA (creación):   Febrero 2021
# FECHA (versión):    23/04/2021 12:43
# CAMBIO MÁS RECIENTE: Cambiado el número de iter y warmup, eliminado (comentado) los gráficos. Incluida una tabla separada en la que se guarden los análisis de convergencia. Arreglado el contador. Solucionado el problema del contador, de nuevo.


############################################
############################################

rm(list=ls())

###################### LIBRERÍAS NECESARIAS  ######################
##                                                               ##
###################################################################
library(MASS) # para mvrnorm
library(lme4) # Linear Mixed-Effects Models
library(nlme) # Linear and Nonlinear Mixed Effects Models (gls)
library(pbkrtest) # Parametric Bootstrap and Kenward Roger Based Methods for Mixed Models Comparison (métodos de estimación con los grados de libertad)
library(lmerTest) # Tests in Linear Mixed Effects Models - para la significación Satterthwaite
library(rstan)
library(brms) # Bayesian Regression Models using 'Stan'
library(bayesplot) # Convergence Diagnostics & Plots
library(ggplot2) # Gráficos
library(tictoc) #to measure simulation time
###################### DIRECTORIO DE TRABAJO  ######################
##                                                                ##
####################################################################

#setwd("/lustre/home/nucifera/CCC/")
setwd("C:/Users/Cris/Documents/Workspace/phdthesis/R")


# -------------------------------- & ---------------------------------- #



####################### DATOS SIMULADOS #####################
##                                                         ##
#############################################################

nrep <- 500 # nº de réplicas
resultados <- matrix(data=NA, nrow = nrep+1, ncol = 316) # Volcamos los resultados en una matriz que tenga el número de filas igual al de las réplicas + 1 y con tantas columnas como tengamos que guardar datos de los parámetros que nos interesen.
resultados_convergencia <- matrix(data=NA, nrow = nrep+1, ncol = 102) #para guardar el análisis de convergencia de las cadenas de Markov Monte Carlo

#Falta por replicar: modelo minimal; te = 2 (alto); MR = 0 y N = 1

NSUJ <- c(5)
MR <- c(5)
#nsuj <- 5 # nº de sujetos del diseño N=1 AB: 3, 5 y 7
#mr <- 20 # nº de medidas repetidas por fase (MR total es mr*2): 10, 20, 30 y 40
#ntotal = nsuj*mr*2 # nº total de observaicones del diseño



##### EFECTOS ALEATORIOS ##
##                       ##
###########################

rho <- 0 # correlación entre las pendientes y las intersecciones (¿otros valores?); en los parciales tiene que ser 0, y en el maximal hay que pensarlo.

varintpend <- 2

SIGMA <- matrix(NA, nrow = 2, ncol = 2) # matriz de varianzas-covarianzas
SIGMA[1,1] <- 2 # varianza de las intersecciones, lo cambiamos a 0, 2 (Moeyaert et al., 2017)
SIGMA[2,2] <- 2 # varianza de las pendientes, lo cambiamos a 0, 2 (Moeyaert et al., 2017)
SIGMA[1,2] <- sqrt(SIGMA[1,1])*sqrt(SIGMA[2,2])*rho #covarianza entre pendientes e intersecciones
SIGMA[2,1] <- SIGMA[1,2] #Es una matriz cuadrada simétrica

MU <- matrix(0, nrow = 1, ncol = 2) #cómo se desvía por término medio un sujeto de gamma00 y gamma10
sigma1 <- 1 # varianza residual, fijada a 1.


##### EFECTOS FIJOS ####
##                    ##
########################
gamma00 <- 5 #línea base fijada a 5
gamma10 <- 2.70 # fijamos el efecto de la intervención bajo la d de Cohen (diferencia individual estandarizada) y los criterios de Ferguson: 0 = no efecto; 1.15 = moderado; 2.70 = grande.


####### BUCLE DE SIMULACIÓN #################
##                                         ##
#############################################

# --------------- Condiciones
tic()
for (v in 1:1){
  nsuj = NSUJ[v] # cambio por cada número de sujetos establecido
  for(w in 1:1){
    mr = MR[w] # cambio para cada número de MR establecidas
    ntotal = mr*nsuj*2 # nº de observaciones totales del diseño

    x <- c() # Variable independiente del estudio con 0s y 1s. (0 = línea base; 1 = tratamiento)
    for (j in 1:nsuj){
      for (i in 1:(2*mr)){
        if (i <= mr){x[i+(j-1)*2*mr]=0}
        if(i > mr){x[i+(j-1)*2*mr]=1}
      }
    }

    #----------- Generación de datos

    for (k in 1:nrep){
      ui = rep(NA, nsuj) # almacenar las intersecciones de los sujetos en la línea base
      us = rep(NA, nsuj) # almacenar las pendientes de los sujetos
      y = rep(NA, ntotal) # VD nivel 1
      id = rep(NA, ntotal) # id del sujeto


      U = matrix(NA, nrow=1, ncol=2)

      for (j in 1:nsuj){
        #ui[j]=rnorm(1,0,sqrt(5))
        #us[j]=rnorm(1,0,sqrt(1))
        #us[j]=0
        U = mvrnorm(1,MU,SIGMA)
        ui[j]=U[1]
        us[j]=U[2]
        for (i in 1:(2*mr)){
          y[i+(j-1)*2*mr] = gamma00 + gamma10*x[i]+rnorm(1, 0, sigma1) #minimal
          id[i+(j-1)*2*mr] = j #identifica la unidad de nivel 2 para los análisis posteriores
        }
      }

      #---------------------- Análisis de los datos

      #  ------------- PERSPECTIVA FRECUENTISTA

      #  - modellmer4: maximal (intersecciones y pendientes aleatorias)
      #  - modellmer3: parcial_2 (intersecciones aleatorias, pendientes fijadas a 0)
      #  - modellmer2: parcial_1 (intersecciones fijadas a 0, pendientes aleatorias)
      #  - modellmer1: minimal (ain pendientes ni intersecciones aleatorias)
      #

      if (k > 1){
        modellmer4 = lmer(y ~ x + (x|id),REML=TRUE) #intesecciones y pendientes, maximal
        modellmer3 = lmer(y ~ x + (1|id),REML=TRUE) #solo intersecciones, parcial_2
        modellmer2 = lmer(y ~ x + (x-1|id), REML = TRUE) #sólo pendientes
        modellmer1 = gls(y ~ x, method = "REML")

        pos = 1;
        pos1 = 1;
        # ----- Extracción de los parámetros de los modelos

        #  MAXIMAL: pendientes e intersecciones aleatorias

        resultados[k,pos]=summary(modellmer4)$coef[1,1];pos=pos+1; #est. g00 from lmer
        resultados[k,pos]=summary(modellmer4)$coef[1,2];pos=pos+1; #se de g00
        resultados[k,pos]=summary(modellmer4)$coef[2,1];pos=pos+1; #est. g10 from lmer
        resultados[k,pos]=summary(modellmer4)$coef[2,2];pos=pos+1; #se de g10
        #resultados[k,pos]=summary(modellmer4)$coef[1,5];pos=pos+1; #valor-p de g00
        resultados[k,pos]=summary(modellmer4)$coef[2,5];pos=pos+1; #valor-p de g10

        resultados[k,pos]=as.data.frame(VarCorr(modellmer4))[4,4];pos=pos+1; #level-1 residual variance from lmer
        resultados[k,pos]=as.data.frame(VarCorr(modellmer4))[1,4];pos=pos+1; #level-2 intercept variance from lmer
        resultados[k,pos]=as.data.frame(VarCorr(modellmer4))[3,4];pos=pos+1; #level-2 int-slope covariance from lmer
        resultados[k,pos]=as.data.frame(VarCorr(modellmer4))[2,4];pos=pos+1; #level-2 slope variance from lmer
        resultados[k,pos] = as.data.frame(VarCorr(modellmer4))[3,5];pos=pos+1; #cor int-slope

        # Índices de ajuste: AIC y BIC (maximal)

        resultados[k,pos]=AIC(modellmer4); pos=pos+1; # AIC MODELO pendientes e intersecciones
        resultados[k,pos]=BIC(modellmer4);pos=pos+1; # BIC MODELO pendientes e intersecciones


        # PARCIAL_2: intersecciones aleatorias

        resultados[k,pos]=summary(modellmer3)$coef[1,1];pos=pos+1; #est. g00 from lmer
        resultados[k,pos]=summary(modellmer3)$coef[1,2];pos=pos+1; #se de g00
        resultados[k,pos]=summary(modellmer3)$coef[2,1];pos=pos+1; #est. g10 from lmer
        resultados[k,pos]=summary(modellmer3)$coef[2,2];pos=pos+1; #se de g10
        #resultados[k,pos]=summary(modellmer3)$coef[1,5];pos=pos+1; #valor-p de g00
        resultados[k,pos]=summary(modellmer3)$coef[2,5];pos=pos+1; #valor-p de g10

        resultados[k,pos]=as.data.frame(VarCorr(modellmer3))[2,4];pos=pos+1; #level-1 residual variance from lmer
        resultados[k,pos]=as.data.frame(VarCorr(modellmer3))[1,4];pos=pos+1; #level-2 intercept variance from lmer

        # Índices de ajuste: AIC y BIC (parcial_2, sólo intersecciones)

        resultados[k,pos]=AIC(modellmer3); pos=pos+1; # AIC MODELO intersecciones
        resultados[k,pos]=BIC(modellmer3);pos=pos+1; # BIC MODELO intersecciones



        # PARCIAL_1: pendientes aleatorias

        resultados[k,pos]=summary(modellmer2)$coef[1,1];pos=pos+1; #est. g00 from lmer
        resultados[k,pos]=summary(modellmer2)$coef[1,2];pos=pos+1; #se de g00
        resultados[k,pos]=summary(modellmer2)$coef[2,1];pos=pos+1; #est. g10 from lmer
        resultados[k,pos]=summary(modellmer2)$coef[2,2];pos=pos+1; #se de g10
        #resultados[k,pos]=summary(modellmer2)$coef[1,5];pos=pos+1; #valor-p de g00
        resultados[k,pos]=summary(modellmer2)$coef[2,5];pos=pos+1; #valor-p de g10

        resultados[k,pos]=as.data.frame(VarCorr(modellmer2))[2,4];pos=pos+1; #level-1 residual variance from lmer
        resultados[k,pos]=as.data.frame(VarCorr(modellmer2))[1,4];pos=pos+1; #level-2 slope variance from lmer

        #Índices de ajuste: AIC y BIC (sólo pendientes)

        resultados[k,pos]=AIC(modellmer2); pos=pos+1 #AIC MODELO PENDIENTES ALEATORIAS
        resultados[k,pos]=BIC(modellmer2); pos=pos+1 #BIC MODELO PENDIENTES ALEATORIAS



        # MINIMAL: sin pendientes ni intersecciones aleatorias
        # Se hace necesario acudir al objeto 'tTable' que devuelve summary(gls)

        resultados[k,pos]=summary(modellmer1)$coef[1];pos=pos+1; #est g00 from gls
        resultados[k, pos]=summary(modellmer1)$tTable[1,2];pos=pos+1; #se of g00
        resultados[k,pos]=summary(modellmer1)$coef[2];pos=pos+1 #est g10 from lmer
        resultados[k,pos]=summary(modellmer1)$tTable[2,2]; pos=pos+1 #se of g10
        #resultados[k,pos]=summary(modellmer1)$tTable[1,4]; pos=pos+1 # valor-p de g00
        resultados[k,pos]=summary(modellmer1)$tTable[2,4]; pos=pos+1 # valor-p de g10

        resultados[k,pos]=summary(modellmer1)$sigma; pos=pos+1 #varianza residual

        #Índices de ajuste: AIC y BIC (ni pendientes ni intersecciones)
        resultados[k,pos]=summary(modellmer1)$AIC;pos=pos+1; # AIC modelo minimal
        resultados[k,pos]=summary(modellmer1)$BIC;pos=pos+1 # BIC modelo minimal


      }



      # ------------- PERSPECTIVA BAYESIANA
      databrms <- cbind(y, x, id)

      #  MAXIMAL: pendientes e intersecciones aleatorias (dos efectos aleatorios)
      #   --- Modelo cauchy 1: half-cauchy(0,10)

      if (k==1){
        fit_mxc1 <- brm(y ~ x + (1+x|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("cauchy(0,10)", class = "sd"),
                                  set_prior("lkj(2)", class="cor"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400,
                        iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_c1 = update(fit_mxc1,newdata=databrms, cores = 4)

        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_c1)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_c1)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_c1)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_c1)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_c1)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_c1)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_c1)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_c1)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_c1)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_c1)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_c1)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # ---- Traceplot para el análisis de convergencia

        #  png("Maximal_cauchy_1.png", height = 500 , width = 700, res = 71)
        #
        #  color_scheme_set("mix-blue-pink")
        #  p <- mcmc_trace(fit_max_c1,  pars = c("b_Intercept", "b_x",
        #                                        "sd_id__Intercept", "sd_id__x",
        #                                        "cor_id__Intercept__x", "sigma"),
        #                  n_warmup = 300,
        #                  facet_args = list(nrow = 2, labeller = label_parsed))
        #  p + facet_text(size = 10)
        #
        # dev.off()


      }


      #   --- Modelo cauchy 2: half-cauchy(0,20)

      if (k==1){
        fit_mxc2 <- brm(y ~ x + (1+x|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("cauchy(0,20)", class = "sd"),
                                  set_prior("lkj(2)", class="cor"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400,
                        iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_c2 = update(fit_mxc2,newdata=databrms, cores = 4)
        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_c2)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_c2)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_c2)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_c2)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_c2)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_c2)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_c2)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_c2)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_c2)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_c2)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_c2)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)


        # ---- Traceplot para el análisis de convergencia

        # png("Maximal_cauchy_2.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_max_c2,  pars = c("b_Intercept", "b_x",
        #                                       "sd_id__Intercept", "sd_id__x",
        #                                       "cor_id__Intercept__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
        #
      }

      #   --- Modelo cauchy 3: half-cauchy(0,50)

      if (k==1){
        fit_mxc3 <- brm(y ~ x + (1+x|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("cauchy(0,50)", class = "sd"),
                                  set_prior("lkj(2)", class="cor"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400,
                        iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_c3 = update(fit_mxc3,newdata=databrms, cores = 4)
        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_c3)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_c3)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_c3)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_c3)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_c3)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_c3)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_c3)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_c3)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_c3)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_c3)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_c3)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # ---- Traceplot para el análisis de convergencia

        # png("Maximal_cauchy_3.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_max_c3,  pars = c("b_Intercept", "b_x", "sd_id__Intercept",
        #                                       "sd_id__x", "cor_id__Intercept__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      #   --- Modelo half-normal 1: half-normal(0,10)

      if (k==1){
        fit_mxn1 <- brm(y ~ x + (1+x|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("normal(0,10)", class = "sd"),
                                  set_prior("lkj(2)", class="cor"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400,
                        iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_n1 = update(fit_mxn1,newdata=databrms, cores = 4)
        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_n1)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_n1)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_n1)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_n1)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_n1)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_n1)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_n1)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_n1)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_n1)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_n1)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_n1)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # ---- Traceplot para el análisis de convergencia

        # png("Maximal_normal_1.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_max_n1,  pars = c("b_Intercept", "b_x", "sd_id__Intercept",
        #                                       "sd_id__x", "cor_id__Intercept__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }


      #   --- Modelo half-normal 2: half-normal(0,20)

      if (k==1){
        fit_mxn2 <- brm(y ~ x + (1+x|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("normal(0,20)", class = "sd"),
                                  set_prior("lkj(2)", class="cor"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400,
                        iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_n2 = update(fit_mxn2,newdata=databrms, cores = 4)
        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_n2)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_n2)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_n2)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_n2)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_n2)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_n2)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_n2)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_n2)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_n2)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_n2)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_n2)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        #
        # png("Maximal_normal_2.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_max_n2,  pars = c("b_Intercept", "b_x", "sd_id__Intercept",
        #                                       "sd_id__x", "cor_id__Intercept__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }


      #   --- Modelo half-normal 3: half-normal(0,50)

      if (k==1){
        fit_mxn3 <- brm(y ~ x + (1+x|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("normal(0,50)", class = "sd"),
                                  set_prior("lkj(2)", class="cor"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400,
                        iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_n3 = update(fit_mxn3,newdata=databrms, cores = 4)
        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_n3)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_n3)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_n3)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_n3)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_n3)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO


        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_n3)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_n3)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_n3)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_n3)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_n3)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_n3)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        #
        # png("Maximal_normal_3.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_max_n3,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sd_id__x",
        #                                       "cor_id__Intercept__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }


      #   --- Modelo uniform: uniform ~(0, 100)

      if (k==1){
        fit_mxuni <- brm(y ~ x + (1+x|id), data = databrms,
                         prior = c(set_prior("normal(0,1000000)", class = "b"),
                                   set_prior("uniform(0,100)", class = "sd"),
                                   set_prior("lkj(2)", class="cor"),
                                   set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                         warmup = 400,
                         iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)
      }

      if (k > 1) {
        fit_max_uni = update(fit_mxuni,newdata=databrms, cores = 4)
        # --- Resultados
        resultados[k,pos]=posterior_summary(fit_max_uni)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_max_uni)[6,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[4,1]^2;pos=pos+1; #est. level-2 slope variance from brms
        resultados[k,pos]=posterior_summary(fit_max_uni)[5,1];pos=pos+1; #est. cor slopes-intercepts from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_max_uni)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_max_uni)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_max_uni)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_max_uni)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_max_uni)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_max_uni)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_max_uni)[5]; pos1=pos1+1; # Gelman-Rubin Statistic for cor slopes-intercepts
        resultados_convergencia[k,pos1] = rhat(fit_max_uni)[6]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Maximal_uniform.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_max_uni,  pars = c("b_Intercept",
        #                                        "b_x", "sd_id__Intercept",
        #                                        "sd_id__x", "cor_id__Intercept__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()

      }



      # PARCIAL 2: sólo intersecciones aleatorias

      # --- Modelo sólo intersecciones half-cauchy 1: half-cauchy(0, 10)
      if (k==1){
        fitp2c1 <- brm(y ~ x + (1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("cauchy(0,10)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_c1 <- update(fitp2c1, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_c1)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_c1)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_c1)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_c1)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_c1)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_c1)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_c1)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_c1)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_c1)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_cauchy_1.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_c1,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo sólo intersecciones half-cauchy 2: half-cauchy(0, 20)
      if (k==1){
        fitp2c2 <- brm(y ~ x + (1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("cauchy(0,20)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_c2 <- update(fitp2c2, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_c2)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_c2)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_c2)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_c2)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_c2)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_c2)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_c2)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_c2)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_c2)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_cauchy_2.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_c2,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()

      }

      # --- Modelo sólo intersecciones half-cauchy 3: half-cauchy(0, 50)
      if (k==1){
        fitp2c3 <- brm(y ~ x + (1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("cauchy(0,50)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_c3 <- update(fitp2c3, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_c3)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_c3)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_c3)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_c3)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_c3)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_c3)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_c3)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_c3)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_c3)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_cauchy_3.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_c3,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }


      # --- Modelo sólo intersecciones half-normal 1: half-normal(0, 10)
      if (k==1){
        fitp2n1 <- brm(y ~ x + (1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("normal(0,10)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_n1 <- update(fitp2n1, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_n1)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_n1)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_n1)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_n1)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_n1)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO


        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_n1)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_n1)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_n1)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_n1)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)


        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_normal_1.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_n1,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo sólo intersecciones half-normal 2: half-normal(0, 20)
      if (k==1){
        fitp2n2 <- brm(y ~ x + (1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("normal(0,20)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_n2 <- update(fitp2n2, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_n2)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_n2)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_n2)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_n2)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_n2)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_n2)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_n2)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_n2)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_n2)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)


        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_normal_2.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_n2,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo sólo intersecciones half-normal 3: half-normal(0, 50)
      if (k==1){
        fitp2n3 <- brm(y ~ x + (1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("normal(0,50)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_n3 <- update(fitp2n3, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_n3)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_n3)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_n3)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_n3)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_n3)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_n3)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_n3)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_n3)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_n3)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_normal_3.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_n3,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo sólo intersecciones uniforme: uniform ~ (0,100)
      if (k==1){
        fitp2uni <- brm(y ~ x + (1|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("uniform(0,100)", class = "sd"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part2_uni <- update(fitp2uni, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part2_uni)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part2_uni)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part2_uni)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part2_uni)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part2_uni)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part2_uni)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part2_uni)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part2_uni)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 intercept variance
        resultados_convergencia[k,pos1] = rhat(fit_part2_uni)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)


        # # ---- Traceplot para el análisis de convergencia
        # png("Partial2_uniforme.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part2_uni,  pars = c("b_Intercept", "b_x", "sd_id__Intercept", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # PARCIAL 1: sólo pendientes aleatorias
      # --- Modelo pendientes aleatorias half-cauchy 1: half-cauchy(0, 10)
      if (k==1){
        fitp1c1 <- brm(y ~ x + (x-1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("cauchy(0,10)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_c1 <- update(fitp1c1, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_c1)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_c1)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_c1)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_c1)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_c1)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_c1)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_c1)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_c1)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_c1)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_cauchy_1.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_c1,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()

      }

      # --- Modelo pendientes aleatorias half-cauchy 2: half-cauchy(0, 20)
      if (k==1){
        fitp1c2 <- brm(y ~ x + (x-1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("cauchy(0,20)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_c2 <- update(fitp1c2, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_c2)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_c2)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_c2)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_c2)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_c2)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO


        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_c2)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_c2)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_c2)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_c2)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_cauchy_2.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_c2,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }


      # --- Modelo pendientes aleatorias half-cauchy 3: half-cauchy(0, 50)
      if (k==1){
        fitp1c3 <- brm(y ~ x + (x-1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("cauchy(0,50)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_c3 <- update(fitp1c3, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_c3)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_c3)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_c3)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_c3)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_c3)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_c3)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_c3)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_c3)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_c3)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_cauchy_3.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_c3,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo pendientes aleatorias half-normal 1: half-normal(0, 10)
      if (k==1){
        fitp1n1 <- brm(y ~ x + (x-1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("normal(0,10)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_n1 <- update(fitp1n1, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_n1)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_n1)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_n1)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_n1)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_n1)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_n1)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_n1)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_n1)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_n1)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_normal_1.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_n1,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo pendientes aleatorias half-normal 2: half-normal(0, 20)
      if (k==1){
        fitp1n2 <- brm(y ~ x + (x-1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("normal(0,20)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_n2 <- update(fitp1n2, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_n2)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_n2)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_n2)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_n2)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_n2)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_n2)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_n2)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_n2)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_n2)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_normal_2.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_n2,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()

      }


      # --- Modelo pendientes aleatorias half-normal 3: half-normal(0, 50)
      if (k==1){
        fitp1n3 <- brm(y ~ x + (x-1|id), data = databrms,
                       prior = c(set_prior("normal(0,1000000)", class = "b"),
                                 set_prior("normal(0,50)", class = "sd"),
                                 set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                       warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_n3 <- update(fitp1n3, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_n3)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_n3)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_n3)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_n3)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_n3)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_n3)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_n3)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_n3)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_n3)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_normal_3.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_n3,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()
      }

      # --- Modelo pendientes aleatorias uniform: uniform(0, 100)
      if (k==1){
        fitp1uni <- brm(y ~ x + (x-1|id), data = databrms,
                        prior = c(set_prior("normal(0,1000000)", class = "b"),
                                  set_prior("uniform(0,100)", class = "sd"),
                                  set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                        warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_part1_uni <- update(fitp1uni, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_part1_uni)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_part1_uni)[4,1]^2;pos=pos+1; #est. de sigma from brms
        resultados[k,pos]=posterior_summary(fit_part1_uni)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms


        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_part1_uni)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_part1_uni)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_part1_uni)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_part1_uni)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_part1_uni)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for level-2 slope variance
        resultados_convergencia[k,pos1] = rhat(fit_part1_uni)[4]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)

        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_uniform.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_uni,  pars = c("b_Intercept", "b_x", "sd_id__x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()

      }


      # MINIMAL: sin efectos aleatorios
      # --- Modelo minimal - al no haber efectos aleatorios, sólo hay un modelo
      if (k==1){
        fitmin <- brm(y ~ x, data = databrms,
                      prior = c(set_prior("normal(0,1000000)", class = "b"),
                                set_prior("inv_gamma(0.001,0.001)", class = "sigma")),
                      warmup = 400, iter = 1000, chains = 2, control = list(adapt_delta = 0.95), cores = 4)

      }
      if (k > 1){
        fit_min_unif <- update(fitmin, newdata = databrms, cores = 4)

        # --- Resultados

        resultados[k,pos]=posterior_summary(fit_min_unif)[1,1];pos=pos+1; #est. g00 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[1,2];pos=pos+1; #se de g00 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[1,3];pos=pos+1; #li de g00 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[1,4];pos=pos+1; #ls de g00 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[2,1];pos=pos+1; #est. g10 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[2,2];pos=pos+1; #se de g10 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[2,3];pos=pos+1; #li de g10 from brms
        resultados[k,pos]=posterior_summary(fit_min_unif)[2,4];pos=pos+1; #ls de g10 from brms

        resultados[k,pos]=posterior_summary(fit_min_unif)[3,1]^2;pos=pos+1; #est. de sigma from brms
        #resultados[k,pos]=posterior_summary(fit_min_unif)[3,1]^2;pos=pos+1; #est. level-2 intercept variance from brms

        # --- Índices de ajuste bayesianos
        waic_v = WAIC(fit_min_unif)
        resultados[k,pos]=waic_v$estimates[3,1];pos=pos+1; #est. WAIC
        loo_v = LOO(fit_min_unif)
        resultados[k,pos]=loo_v$estimates[3,1];pos=pos+1; #est. LOO

        # --- Análisis de convergencia
        resultados_convergencia[k,pos1] = rhat(fit_min_unif)[1]; pos1=pos1+1; # Gelman-Rubin Statistic for g00
        resultados_convergencia[k,pos1] = rhat(fit_min_unif)[2]; pos1=pos1+1; # Gelman-Rubin Statistic for g10
        resultados_convergencia[k,pos1] = rhat(fit_min_unif)[3]; pos1=pos1+1; # Gelman-Rubin Statistic for sigma (level-1 residual variance)


        # # ---- Traceplot para el análisis de convergencia
        # png("Partial1_normal_1.png", height = 500 , width = 700, res = 71)
        #
        # color_scheme_set("mix-blue-pink")
        # p <- mcmc_trace(fit_part1_n1,  pars = c("b_Intercept", "b_x", "sigma"),
        #                 n_warmup = 300,
        #                 facet_args = list(nrow = 2, labeller = label_parsed))
        # p + facet_text(size = 10)
        #
        # dev.off()


        print(k)
        print(pos)
        print(pos1)

      }
    }

    mrs = toString(mr);nsujs = toString(nsuj); g10s = toString(gamma10); varintpendi = toString(varintpend)
    filenombre = paste("MINIMAL", g10s, sep = " ")
    filenombre = paste(filenombre, mrs,sep = " ")
    filenombre = paste(filenombre,nsujs,sep = " ")
    filenombre = paste(filenombre, varintpendi, sep = " ")
    filenombre = paste(filenombre,".dat",sep="")

    write.table(resultados, filenombre,sep = ";",row.names = FALSE)


    filenombre = paste("MINIMAL convergencia", g10s, sep = " ")
    filenombre = paste(filenombre, mrs,sep = " ")
    filenombre = paste(filenombre,nsujs,sep = " ")
    filenombre = paste(filenombre, varintpendi, sep = " ")
    filenombre = paste(filenombre,".dat",sep="")
    write.table(resultados_convergencia, filenombre,sep = ";",row.names = FALSE)

  }
}
toc()
