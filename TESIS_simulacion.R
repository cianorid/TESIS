###############################################################################
# Ejemplo completo de MCMC para un modelo SV
#
#   y_t = e_t,    e_t ~ N(0, exp(h_t))
#   h_t = mu + phi_h * (h_{t-1} - mu) + sigma_h * eta_t,   eta_t ~ N(0,1)
#
###############################################################################

# Cargar librerías necesarias
library(MASS)    # Para mvrnorm
library(coda)    # Para diagnósticos MCMC
library(ggplot2)

source("TESIS_funciones.R")

# ---------------------------------------------------------------------
# 1) Función para simular datos SV
# ---------------------------------------------------------------------
simular_SV <- function(T, mu, phi_h, sigma_h) {
  y <- numeric(T)
  h <- numeric(T)
  y[1] <- 0      # Inicializar y_1
  h[1] <- mu     # Inicializar h_1
  
  for (t in 2:T) {
    h[t] <- mu + phi_h * (h[t-1] - mu) + sigma_h * rnorm(1)
    y[t] <- rnorm(1, mean = 0, sd = sqrt(exp(h[t])))
  }
  
  return(list(y = y, h = h))
}

# ---------------------------------------------------------------------
# 2) Ejemplo de uso
# ---------------------------------------------------------------------
set.seed(124)
# Parámetros verdaderos del modelo SV
mu_true <- 3
phi_h_true <- 0.5
sigma_h_true <- 0.5
N_mcmc <- 50000

# Simular datos
sim <- simular_SV(T = 200, mu = mu_true, phi_h = phi_h_true, sigma_h = sigma_h_true)
y_data <- sim$y
h_true <- sim$h  # Volatilidades reales

# Visualizar los datos simulados
df <- data.frame(Tiempo = 1:length(y_data), y_t = y_data)
ggplot(df, aes(x = Tiempo, y = y_t)) +
  geom_line(color = "black", size = 0.3) +
  theme_minimal(base_size = 14) +
  labs(title = expression("Datos Simulados" ~ y[t]), x = "Tiempo", y = expression(y[t])) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))




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

# Ejecutar el MCMC
res_mcmc <- SV_Gibbs(y = y_data, N_mcmc = N_mcmc,
                     proposal_sd_h = 1,
                     prior_sigma_shape = 2.5, prior_sigma_rate = 0.025,
                     b0 = b0, Sigma0 = Sigma0,
                     mu_init = 0, phi_h_init = 0.5, sigma_h_init = 0.1)

cat("Tasa de aceptación conjunta (mu, phi_h):", res_mcmc$tasa_acept_mu_phi, "\n")
cat("Tasa de aceptación (h_t):", res_mcmc$tasa_acept_h, "\n")


library(ggplot2)
library(gridExtra)
library(coda)





# ---------------------------------------------------------------------
# Visualización de Trazas de MCMC (Antes de Thinning y Burn-in)
# ---------------------------------------------------------------------
plot_trace_pre_thinning <- function(chain, param_name) {
  ggplot(data.frame(Iteración = 1:length(chain), Valor = chain), aes(x = Iteración, y = Valor)) +
    geom_line(color = "darkorange") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Traza MCMC para", param_name, "(Antes de Thinning)"), 
         y = param_name, x = "Iteración")
}

# plot_mu_pre <- plot_trace_pre_thinning(res_mcmc$mu, expression(mu))
# plot_phi_h_pre <- plot_trace_pre_thinning(res_mcmc$phi_h, expression(phi[h]))
# plot_sigma_h_pre <- plot_trace_pre_thinning(res_mcmc$sigma_h, expression(sigma[h]))

# Mostrar las trazas juntas
# grid.arrange(plot_mu_pre, nrow = 1)
# grid.arrange(plot_phi_h_pre, nrow = 1)
# grid.arrange(plot_sigma_h_pre, nrow = 1)


# ---------------------------------------------------------------------
# 1) Autocorrelación Inicial
# ---------------------------------------------------------------------
# acf(res_mcmc$mu, main = expression("Autocorrelación de " * mu), lag.max = 40)
# acf(res_mcmc$phi_h, main = expression("Autocorrelación de " * phi[h]), lag.max = 40)
# acf(res_mcmc$sigma_h, main = expression("Autocorrelación de " * sigma[h]), lag.max = 40)

# ---------------------------------------------------------------------
# 2) Aplicar Thinning (si es necesario)
# ---------------------------------------------------------------------
thinning_step <-200
burn_in <- 3000
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
# 4) Cálculo de promedios posteriores
# ---------------------------------------------------------------------
mu_posterior_mean <- mean(mu, na.rm = TRUE)
phi_h_posterior_mean <- mean(phi_h, na.rm = TRUE)
sigma_h_posterior_mean <- mean(sigma_h, na.rm = TRUE)

cat("Promedios posteriores:\n")
cat("mu =", round(mu_posterior_mean, 3), "\n")
cat("phi_h =", round(phi_h_posterior_mean, 3), "\n")
cat("sigma_h =", round(sigma_h_posterior_mean, 3), "\n")

# ---------------------------------------------------------------------
# 5) Visualización de Trazas de MCMC
# ---------------------------------------------------------------------
plot_trace <- function(chain, true_value, param_name) {
  ggplot(data.frame(Iteración = 1:length(chain), Valor = chain), aes(x = Iteración, y = Valor)) +
    geom_line(color = "steelblue") +
    geom_hline(yintercept = true_value, color = "red", linetype = "dashed", size = 1) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Traza MCMC para", param_name), y = param_name, x = "Iteración")
}

plot_mu <- plot_trace(mu, mu_true, expression(mu))
plot_phi_h <- plot_trace(phi_h, phi_h_true, expression(phi[h]))
plot_sigma_h <- plot_trace(sigma_h, sigma_h_true, expression(sigma[h]))

grid.arrange(plot_mu, plot_phi_h, plot_sigma_h, nrow = 1)

# ---------------------------------------------------------------------
# 6) Resúmenes Posteriores de h_t
# ---------------------------------------------------------------------
h_posterior_samples <- res_mcmc$h[(burn_in + 1):N_mcmc, ]
h_posterior_mean <- apply(h_posterior_samples, 2, mean)
h_posterior_lower <- apply(h_posterior_samples, 2, quantile, probs = 0.025)
h_posterior_upper <- apply(h_posterior_samples, 2, quantile, probs = 0.975)

tiempo <- 1:length(h_true)
subset_start <- 1
subset_end <- 200
tiempo_subset <- tiempo[subset_start:subset_end]
h_true_subset <- h_true[subset_start:subset_end]
h_mean_subset <- h_posterior_mean[subset_start:subset_end]
h_lower_subset <- h_posterior_lower[subset_start:subset_end]
h_upper_subset <- h_posterior_upper[subset_start:subset_end]

plot_h_comparison <- ggplot() +
  geom_line(aes(x = tiempo_subset, y = h_true_subset), color = "black", size = 0.5) +
  geom_line(aes(x = tiempo_subset, y = h_mean_subset), color = "blue", size = 0.5) +
  geom_ribbon(aes(x = tiempo_subset, ymin = h_lower_subset, ymax = h_upper_subset), fill = "red", alpha = 0.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Volatilidades Reales vs. Estimadas con IC del 95%",
       x = "Tiempo", y = expression(h[t]))

print(plot_h_comparison)

# ---------------------------------------------------------------------
# 7) Histogramas de mu, phi_h y sigma_h
# ---------------------------------------------------------------------
plot_hist <- function(chain, true_value, param_name) {
  ggplot(data.frame(Valor = chain), aes(x = Valor)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = true_value, color = "red", linetype = "dashed", size = 1) +
    theme_minimal(base_size = 14) +
    labs(title = paste( param_name), x = param_name, y = "Frecuencia")
}

plot_hist_mu <- plot_hist(mu, mu_true, expression(mu))
plot_hist_phi_h <- plot_hist(phi_h, phi_h_true, expression(phi[h]))
plot_hist_sigma_h <- plot_hist(sigma_h, sigma_h_true, expression(sigma[h]))

grid.arrange(plot_hist_mu, plot_hist_phi_h, plot_hist_sigma_h, nrow = 1)



