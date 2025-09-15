###############################################################################
# Ejemplo completo de MCMC para un modelo SV
#
#   y_t = e_t,    e_t ~ N(0, exp(h_t))
#   h_t = mu + phi_h * (h_{t-1} - mu) + sigma_h * eta_t,   eta_t ~ N(0,1)
#
###############################################################################

# Cargar librerías necesarias
library(MASS)     # Para mvrnorm
library(coda)     # Para diagnósticos MCMC
library(ggplot2)  # Para graficos
library(quantmod) # Para extraer datos financieros
library(gridExtra)
library(plotly)




source("TESIS_funciones.R")

# ---------------------------------------------------------------------
# 1) Datos verdaderos
# ---------------------------------------------------------------------

# --- Datos de Entrenamiento ---
symbol <- 'BTC-USD'
start_date <- '2022-01-01'
end_date <- '2024-01-01'

# Descargar datos de Yahoo Finance para el período de entrenamiento
datos <- getSymbols(symbol, src = 'yahoo', from = start_date, to = end_date, auto.assign = FALSE)

# Extraer precios de cierre y calcular retornos logarítmicos
y_data <- as.numeric(Cl(datos))
y_data <- diff(log(y_data))


# Visualizar los datos simulados
df <- data.frame(Tiempo = 1:length(y_data), y_t = y_data)
ggplot(df, aes(x = Tiempo, y = y_t)) +
  geom_line(color = "black", size = 0.3) +
  theme_minimal(base_size = 14) +
  labs(title = expression("Retornos" ~ y[t]), x = "Tiempo", y = expression(y[t])) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))



# Calcular media y desviación estándar del vector
media <- mean(y_data, na.rm = TRUE)
sd_val <- sd(y_data, na.rm = TRUE)

# Convertir el vector en un data frame
df <- data.frame(valor = y_data)

# Crear el histograma con la densidad y la curva de la normal superpuesta
ggplot(df, aes(x = valor)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = media, sd = sd_val), color = "red", size = 1) +
  theme_minimal(base_size = 14) +
  labs(title = "Histograma con distribución normal superpuesta",
       x = "Valor", y = "Densidad")



# -------------------------------
# 3) MCMC
# -------------------------------


# Especificaciones para el prior conjunto de (mu, phi_h)
# Supongamos que, a priori, se espera que mu esté cerca de 0 y phi_h cerca de 0.5.
# Dado que en la reparametrización beta1 = mu*(1-phi_h) y beta2 = phi_h,
# definimos:
prior_mu <- 0
prior_phi <- 0.8
b0 <- c(prior_mu * (1 - prior_phi), prior_phi)  # Por ejemplo: c(0, 0.5)
Sigma0 <- diag(c(100, 1))  # Matriz de covarianza (inversa de la matriz de precisión A)
N_mcmc <- 50000
set.seed(1234)

# Ejecutar el MCMC
res_mcmc <- SV_Gibbs(y = y_data, N_mcmc = N_mcmc,
                     proposal_sd_h = 1.5,
                     prior_sigma_shape = 1, prior_sigma_rate = 0.02,
                     b0 = b0, Sigma0 = Sigma0,
                     mu_init = 0, phi_h_init = 0.001, sigma_h_init = 0.1)

cat("Tasa de aceptación conjunta (mu, phi_h):", res_mcmc$tasa_acept_mu_phi, "\n")
cat("Tasa de aceptación (h_t):", res_mcmc$tasa_acept_h, "\n")





# Visualización de Trazas de MCMC (Antes de Thinning y Burn-in)
# plot_trace_pre_thinning <- function(chain, param_name) {
#   ggplot(data.frame(Iteración = 1:length(chain), Valor = chain), aes(x = Iteración, y = Valor)) +
#     geom_line(color = "darkorange") +
#     theme_minimal(base_size = 14) +
#     labs(title = paste("Traza MCMC para", param_name, "(Antes de Thinning)"), 
#          y = param_name, x = "Iteración")
# }
# 
# plot_mu_pre <- plot_trace_pre_thinning(res_mcmc$mu, expression(mu))
# plot_phi_h_pre <- plot_trace_pre_thinning(res_mcmc$phi_h, expression(phi[h]))
# plot_sigma_h_pre <- plot_trace_pre_thinning(res_mcmc$sigma_h, expression(sigma[h]))

# Imprimir cada gráfico por separado
# print(plot_mu_pre)
# print(plot_phi_h_pre)
# print(plot_sigma_h_pre)


# ---------------------------------------------------------------------
# 1) Autocorrelación Inicial
# ---------------------------------------------------------------------
acf(res_mcmc$mu, main = expression("Autocorrelación de " * mu), lag.max = 40)
acf(res_mcmc$phi_h, main = expression("Autocorrelación de " * phi[h]), lag.max = 40)
acf(res_mcmc$sigma_h, main = expression("Autocorrelación de " * sigma[h]), lag.max = 40)

# ---------------------------------------------------------------------
# 2) Aplicar Thinning (si es necesario)
# -------------------------------------------------------------------
thinning_step <- 120
burn_in <- 5000
mu <- res_mcmc$mu[seq((burn_in + 1), N_mcmc, by = thinning_step)]
phi_h <- res_mcmc$phi_h[seq((burn_in + 1), N_mcmc, by = thinning_step)]
sigma_h <- res_mcmc$sigma_h[seq((burn_in + 1), N_mcmc, by = thinning_step)]

# ---------------------------------------------------------------------
# 1) Autocorrelación despues de Thinning
# ---------------------------------------------------------------------
acf(mu, main = expression("Autocorrelación de " * mu), lag.max = 40)
acf(phi_h, main = expression("Autocorrelación de " * phi[h]), lag.max = 40)
acf(sigma_h, main = expression("Autocorrelación de " * sigma[h]), lag.max = 40)

# ---------------------------------------------------------------------
# 3) Calcular el tamaño de muestra efectivo (ESS)
# ---------------------------------------------------------------------
ess_mu <- effectiveSize(mu)
ess_phi_h <- effectiveSize(phi_h)
ess_sigma_h <- effectiveSize(sigma_h)

cat("Effective Sample Size (ESS):\n")
cat("ESS para mu =", round(ess_mu, 2), "\n")
cat("ESS para phi_h =", round(ess_phi_h, 2), "\n")
cat("ESS para sigma_h =", round(ess_sigma_h, 2), "\n")

# ---------------------------------------------------------------------
# 4) Cálculo de promedios posteriores e intervalos de credibilidad
# ---------------------------------------------------------------------

# Promedios posteriores
mu_posterior_mean      <- mean(mu, na.rm = TRUE)
phi_h_posterior_mean   <- mean(phi_h, na.rm = TRUE)
sigma_h_posterior_mean <- mean(sigma_h, na.rm = TRUE)

# Intervalos de credibilidad (95%)
mu_CI      <- quantile(mu,      probs = c(0.025, 0.975), na.rm = TRUE)
phi_h_CI   <- quantile(phi_h,   probs = c(0.025, 0.975), na.rm = TRUE)
sigma_h_CI <- quantile(sigma_h, probs = c(0.025, 0.975), na.rm = TRUE)

# Mostrar resultados
cat("Resultados posteriores:\n")
cat("mu =", round(mu_posterior_mean, 3), 
    "IC95% [", round(mu_CI[1], 3), ",", round(mu_CI[2], 3), "]\n")
cat("phi_h =", round(phi_h_posterior_mean, 3), 
    "IC95% [", round(phi_h_CI[1], 3), ",", round(phi_h_CI[2], 3), "]\n")
cat("sigma_h =", round(sigma_h_posterior_mean, 3), 
    "IC95% [", round(sigma_h_CI[1], 3), ",", round(sigma_h_CI[2], 3), "]\n")


# ---------------------------------------------------------------------
# 5) Visualización de Trazas de MCMC
# ---------------------------------------------------------------------
plot_trace <- function(chain, param_name) {
  ggplot(data.frame(Iteración = 1:length(chain), Valor = chain), aes(x = Iteración, y = Valor)) +
    geom_line(color = "steelblue") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Traza MCMC para", param_name), y = param_name, x = "Iteración")
}

plot_mu <- plot_trace(mu, expression(mu))
plot_phi_h <- plot_trace(phi_h, expression(phi[h]))
plot_sigma_h <- plot_trace(sigma_h, expression(sigma[h]))

grid.arrange(plot_mu, plot_phi_h, plot_sigma_h, nrow = 1)


# ---------------------------------------------------------------------
# 6) Resúmenes Posteriores de h_t
# ---------------------------------------------------------------------
h_posterior_samples <- res_mcmc$h[(burn_in + 1):N_mcmc, ]
h_posterior_mean <- apply(h_posterior_samples, 2, mean)
h_posterior_lower <- apply(h_posterior_samples, 2, quantile, probs = 0.025)
h_posterior_upper <- apply(h_posterior_samples, 2, quantile, probs = 0.975)

tiempo <- 1:length(h_posterior_mean)
subset_start <- 1
subset_end <- length(y_data)
tiempo_subset <- tiempo[subset_start:subset_end]
h_mean_subset <- h_posterior_mean[subset_start:subset_end]
h_lower_subset <- h_posterior_lower[subset_start:subset_end]
h_upper_subset <- h_posterior_upper[subset_start:subset_end]


plot_h_estimation <- ggplot() +
  geom_line(aes(x = tiempo_subset, y = h_mean_subset), color = "blue", size = 0.8) +
  geom_ribbon(aes(x = tiempo_subset, ymin = h_lower_subset, ymax = h_upper_subset), fill = "red", alpha = 0.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Estimación de h con IC del 95%",
       x = "Tiempo", y = expression(h[t]))

print(plot_h_estimation)



# Supongamos que y_data y h_posterior_mean son vectores de la misma longitud
df <- data.frame(
  Tiempo = 1:length(y_data),
  Retornos = y_data,
  Volatilidad = exp(h_posterior_mean)
)

# Gráfico de Retornos
p1 <- ggplot(df, aes(x = Tiempo, y = Retornos)) +
  geom_line(color = "steelblue", size = 0.7) +
  labs(title = "Retornos", x = "Tiempo", y = "Retornos") +
  theme_minimal(base_size = 14)

# Gráfico de Volatilidad
p2 <- ggplot(df, aes(x = Tiempo, y = Volatilidad)) +
  geom_line(color = "red", size = 0.7) +
  labs(title = "Volatilidad Estocástica", x = "Tiempo", y = "Volatilidad (exp)") +
  theme_minimal(base_size = 14)

# Organizar los gráficos en una sola ventana, uno debajo del otro
grid.arrange(p1, p2, nrow = 2)



# ---------------------------------------------------------------------
# 7) Histogramas de mu, phi_h y sigma_h (sin valores verdaderos)
# ---------------------------------------------------------------------
plot_hist <- function(chain, param_name) {
  ggplot(data.frame(Valor = chain), aes(x = Valor)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    theme_minimal(base_size = 14) +
    labs(title = paste(param_name), x = param_name, y = "Frecuencia")
}

plot_hist_mu <- plot_hist(mu, expression(mu))
plot_hist_phi_h <- plot_hist(phi_h, expression(phi[h]))
plot_hist_sigma_h <- plot_hist(sigma_h, expression(sigma[h]))

grid.arrange(plot_hist_mu, plot_hist_phi_h, plot_hist_sigma_h, nrow = 1)




# ----------------------------------------------------------------------------
# 8) Predicion mediante filtro de particulas, con parametros estimados con SV
# ----------------------------------------------------------------------------



# --- Datos para Prediccion ---
symbol <- 'BTC-USD'
start_date <- '2024-01-02'
end_date <- '2025-01-01'

# Descargar datos de Yahoo Finance para el período de entrenamiento
datos <- getSymbols(symbol, src = 'yahoo', from = start_date, to = end_date, auto.assign = FALSE)

# Extraer precios de cierre y calcular retornos logarítmicos
y_data <- as.numeric(Cl(datos))
y_data <- diff(log(y_data))



source("TESIS_filtro_de_particulas.R")


# Aplicar el Filtro de Partículas con los parámetros estimados
N_particulas <- 1000
h_est <- filtro_particulas(y_data, N_particulas, mu_posterior_mean, phi_h_posterior_mean, sigma_h_posterior_mean )

diagnostico_prediccion_chi <- function(y_data, h_estimado, phi_h, mu, sigma_h, N_particulas = 1000) {
  T <- length(y_data)
  
  # Probabilidades empíricas de predicción
  prob_prediccion <- numeric(T-1)
  
  for(t in 1:(T-1)) {
    h_t1_samples <- mu + phi_h * (h_estimado[t] - mu) + sigma_h * rnorm(N_particulas)
    valores_chi <- (y_data[t+1]^2) / exp(h_t1_samples)
    prob_particula <- pchisq(valores_chi, df = 1)
    prob_prediccion[t] <- mean(prob_particula)
  }
  
  # Transformación a normal estándar
  normal_prediccion <- qnorm(prob_prediccion)
  
  # Tests
  box_test <- Box.test(prob_prediccion, type = "Ljung-Box", lag = 10)
  normalidad_test <- shapiro.test(normal_prediccion)
  modelo_hetero <- lm(normal_prediccion^2 ~ seq_along(normal_prediccion))
  r2_hetero <- summary(modelo_hetero)$r.squared
  
  # === Gráficos individuales ===
  
  # Histograma u_t
  hist(prob_prediccion, breaks = 100, probability = TRUE, col = "lightblue",
       main = "Histograma de u_t", xlab = "u_t", ylab = "Densidad")
  curve(dunif(x, 0, 1), add = TRUE, col = "red", lwd = 2)
  
  # Densidad u_t
  plot(density(prob_prediccion), col = "blue", lwd = 2,
       main = "Densidad de u_t", xlab = "u_t", ylab = "Densidad")
  abline(h = 1, col = "red", lty = 2)
  
  # Histograma n_t
  hist(normal_prediccion, breaks = 50, probability = TRUE, col = "lightgreen",
       main = "Histograma de n_t (Transformación Normal)", xlab = "n_t", ylab = "Densidad", ylim = c(0, 0.6))
  curve(dnorm(x, 0, 1), add = TRUE, col = "red", lwd = 2)
  
  # Densidad n_t
  plot(density(normal_prediccion), col = "blue", lwd = 2,
       main = "Densidad de n_t", xlab = "n_t", ylab = "Densidad")
  curve(dnorm(x, 0, 1), add = TRUE, col = "red", lwd = 2)
  
  return(list(
    u_t = prob_prediccion,
    n_t = normal_prediccion,
    Box_Ljung_pvalor = box_test$p.value,
    Normalidad_pvalor = normalidad_test$p.value,
    Heterocedasticidad_R2 = r2_hetero
  ))
}


# Diagnóstico en la predicción
result_prediction_chi <- diagnostico_prediccion_chi(y_data, h_est, phi_h, mu, sigma_h)
print("Diagnóstico en la predicción (usando Chi-cuadrado escalada con Y^2):")
print(result_prediction_chi)



# -------------------------------
# Densiades evolutivas de y_t
# -------------------------------


graficar_densidades_evolutivas_facets <- function(h_estimado,
                                                   tiempos = 15,
                                                   rango_y = c(-0.15, 0.15)) {
   T      <- length(h_estimado)
   idx_t  <- round(seq(1, T, length.out = tiempos))
   y_vals <- seq(rango_y[1], rango_y[2], length.out = 500)

   # construir data.frame largo
   df_list <- lapply(idx_t, function(t) {
     var_t <- exp(h_estimado[t])
     data.frame(t   = factor(t),
                y   = y_vals,
                pdf = dnorm(y_vals, 0, sqrt(var_t)))
   })
   datos <- do.call(rbind, df_list)

   ggplot(datos, aes(x = y, y = pdf)) +
     geom_line(colour = "steelblue", linewidth = 1) +
     facet_wrap(~ t, ncol = 3, scales = "free_y") +
     labs(title = "Evolución de la distribución condicional de y_t",
          x = "y", y = "Densidad") +
     theme_bw()
 }

# uso:
graficar_densidades_evolutivas_facets(h_est)



# ------------------------------------------------
# Densiades evolutivas de y_t juntas en un grafico
# ------------------------------------------------

graficar_densidades_evolutivas_gg <- function(h_estimado,
                                              tiempos = 15,
                                              rango_y = c(-0.15, 0.15)) {
  # --- 1) Seleccionar los índices de tiempo que vamos a mostrar
  T          <- length(h_estimado)
  idx_t      <- round(seq(1, T, length.out = tiempos))   # 15 tiempos equiespaciados
  
  # --- 2) Mallado de valores de y
  y_vals <- seq(rango_y[1], rango_y[2], length.out = 500)
  
  # --- 3) Construir data.frame largo: una fila por (t, y)
  df_list <- lapply(idx_t, function(t) {
    var_t <- exp(h_estimado[t])
    data.frame(
      t     = t,
      y     = y_vals,
      dens  = dnorm(y_vals, 0, sqrt(var_t))
    )
  })
  datos <- do.call(rbind, df_list)
  
  # --- 4) Graficar con ggplot2
  ggplot(datos, aes(x = y, y = dens, colour = t, group = t)) +
    geom_line(size = 1.1) +
    scale_colour_gradient(low = "blue", high = "red") +
    labs(title = "Evolución de la distribución condicional de y_t",
         x = "y", y = "Densidad", colour = "Tiempo (t)") +
    theme_bw()
}

# Ejemplo de uso:
graficar_densidades_evolutivas_gg(h_est)




 # -------------------------------
# Mezcla de las Normales estimads
# -------------------------------


graficar_distribucion_marginal_y <- function(h_est, y_vals = seq(-0.15, 0.15, length.out = 1000)) {
  # Matriz de densidades: cada columna corresponde a un h
  densidades <- sapply(h_est, function(h) dnorm(y_vals, mean = 0, sd = sqrt(exp(h))))
  
  # Promedio por fila (mezcla empírica)
  densidad_marginal <- rowMeans(densidades)
  
  # Data frame para ggplot
  df <- data.frame(y = y_vals, densidad = densidad_marginal)
  
  # Gráfico con ggplot2
  ggplot(df, aes(x = y, y = densidad)) +
    geom_line(color = "darkblue", linewidth = 1.2) +
    labs(title = "Distribución marginal de Y (mezcla de normales)",
         x = "y", y = "Densidad") +
    theme_minimal()
}




graficar_distribucion_marginal_y(h_est)


# ----------------------------------------------------
# Mezcla de las Normales estimadas con los datos reales
# -----------------------------------------------------


comparar_mezcla_con_datos_ggplot <- function(y_data,
                                             h_est,
                                             y_grid = seq(min(y_data), max(y_data), length.out = 1000),
                                             nbins = 50) {
  # 1) Mezcla empírica p(y)
  dens_mat <- sapply(h_est, function(h) dnorm(y_grid, 0, sqrt(exp(h))))
  dens_mix <- rowMeans(dens_mat)  # p_mezcla(y)
  
  # 2) Crear data frame para mezcla
  df_mezcla <- data.frame(y = y_grid, densidad = dens_mix)
  
  # 3) Crear data frame para histograma de los datos reales
  df_hist <- data.frame(y = y_data)
  
  # 4) Graficar con ggplot2
  p <- ggplot() +
    geom_histogram(data = df_hist, aes(x = y, y = ..density..),
                   bins = nbins, fill = "grey70", color = "grey50", alpha = 0.6) +
    geom_line(data = df_mezcla, aes(x = y, y = densidad),
              color = "blue", size = 1.2) +
    labs(title = "Retornos reales vs. mezcla empírica",
         x = "y", y = "Densidad") +
    theme_minimal()
  
  # 5) Test KS opcional
  cdf_mix <- cumsum(dens_mix) / sum(dens_mix)
  cdf_fun <- approxfun(y_grid, cdf_mix, yleft = 0, yright = 1)
  ks_res  <- ks.test(y_data, cdf_fun)
  
  print(p)
  
  invisible(list(dens_mix = dens_mix,
                 ks_pvalor = ks_res$p.value,
                 ks_stat = ks_res$statistic))
}


resultado_ks <- comparar_mezcla_con_datos_ggplot(y_data, h_est)

print(resultado_ks)


