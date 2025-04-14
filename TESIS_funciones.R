
# ---------------------------------------------------------------------
# 1) Funciones para actualizar parámetros
# ---------------------------------------------------------------------

# 1.1) Actualizar h_t (actualización single-site Metropolis)
actualizar_h_t <- function(t, h, y, mu, phi_h, sigma_h, proposal_sd) {
  # Proponer nuevo h_t
  h_propuesto <- rnorm(1, mean = h[t], sd = proposal_sd)
  
  # Función para calcular el log-posterior de h_t
  logpost_h_t <- function(ht) {
    T <- length(y)
    ll <- 0
    # Contribución de la verosimilitud para y_t
    ll <- ll - 0.5 * ((y[t]^2)/exp(ht) + ht)
    
    # Contribución del proceso AR(1)
    if (t > 1) {
      media_ht <- mu + phi_h * (h[t-1] - mu)
      ll <- ll - 0.5*((ht - media_ht)^2)/sigma_h^2 - 0.5*log(sigma_h^2)
    } else {
      media_ht <- mu
      ll <- ll - 0.5*((ht - media_ht)^2)/sigma_h^2 - 0.5*log(sigma_h^2)
    }
    if (t < T) {
      media_ht1 <- mu + phi_h*(ht - mu)
      ll <- ll - 0.5*((h[t+1]-media_ht1)^2)/sigma_h^2 - 0.5*log(sigma_h^2)
    }
    return(ll)
  }
  
  logpost_actual <- logpost_h_t(h[t])
  logpost_propuesto <- logpost_h_t(h_propuesto)
  log_ratio <- logpost_propuesto - logpost_actual
  
  if (log(runif(1)) < log_ratio) {
    return(list(h_t_nuevo = h_propuesto, aceptado = 1))
  } else {
    return(list(h_t_nuevo = h[t], aceptado = 0))
  }
}

# 1.2) Actualizar conjuntamente (mu, phi_h) mediante Gibbs
# Se reparametriza: beta1 = mu*(1-phi_h) y beta2 = phi_h.
# Para t = 2,...,T se tiene: h_t = beta1 + beta2 * h_{t-1} + error, error ~ N(0, sigma_h^2).
# Se asume un prior: beta ~ N(b0, sigma_h^2 * Sigma0)
actualizar_mu_phi_joint <- function(h, sigma_h, b0, Sigma0) {
  T <- length(h)
  # Construir la matriz de diseño y vector respuesta para t = 2,...,T
  X <- cbind(rep(1, T-1), h[1:(T-1)])
  y_reg <- h[2:T]
  
  # Posterior conjugada: 
  # Varianza posterior: sigma_h^2 * (X'X + Sigma0^{-1})^{-1}
  Sigma_post_inv <- t(X) %*% X + solve(Sigma0)
  Sigma_post <- solve(Sigma_post_inv)
  m_beta <- Sigma_post %*% ( t(X) %*% y_reg + solve(Sigma0) %*% b0 )
  
  # Extraer beta = (beta1, beta2)
  beta <- mvrnorm(1, mu = as.vector(m_beta), Sigma = sigma_h^2 * Sigma_post)
  
  phi_new <- beta[2]
  # Para recuperar mu se usa: beta1 = mu*(1-phi_new)  => mu = beta1/(1-phi_new)
  if (abs(1 - phi_new) < 1e-6) {
    mu_new <- b0[1]  # en caso de inestabilidad, se usa el valor del prior
  } else {
    mu_new <- beta[1] / (1 - phi_new)
  }
  return(c(mu_new, phi_new))
}

# 1.3) Actualizar sigma_h (con posterior Inverse-Gamma)
actualizar_sigma_h <- function(h, mu, phi_h, prior_sigma_shape, prior_sigma_rate) {
  T <- length(h)
  residuals <- h[2:T] - mu - phi_h * (h[1:(T-1)] - mu)
  RSS <- sum(residuals^2)
  shape_post <- prior_sigma_shape + (T - 1) / 2
  rate_post <- prior_sigma_rate + RSS / 2
  sigma2_nuevo <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  sigma_h_nuevo <- sqrt(sigma2_nuevo)
  return(sigma_h_nuevo)
}

# ---------------------------------------------------------------------
# 2) Función principal de Gibbs Sampling
# ---------------------------------------------------------------------
# Se actualizan:
#   1. Los volátiles h[t] (cada uno con MH)
#   2. Conjuntamente (mu,phi_h) mediante Gibbs, usando la reparametrización:
#         beta1 = mu*(1-phi_h)   y   beta2 = phi_h.
#       Se requiere especificar el vector y matriz de covarianza del prior para beta.
#   3. sigma_h se actualiza mediante un paso Inverse-Gamma.
SV_Gibbs <- function(y, N_mcmc, proposal_sd_h,
                     prior_sigma_shape, prior_sigma_rate,
                     b0, Sigma0,   # Prior para beta = (mu*(1-phi_h), phi_h)
                     mu_init, phi_h_init, sigma_h_init) {
  T <- length(y)
  # Contenedores para guardar las muestras
  mu_save <- numeric(N_mcmc)
  phi_h_save <- numeric(N_mcmc)
  sigma_h_save <- numeric(N_mcmc)
  h_save <- matrix(0, nrow = N_mcmc, ncol = T)
  
  # Inicialización
  mu <- mu_init
  phi_h <- phi_h_init
  sigma_h <- sigma_h_init
  h <- rep(log(var(y + 1e-6)), T)  # Inicialización de h_t (se puede ajustar)
  
  acept_h_total <- 0
  acept_mu_phi <- 0
  
  for (iter in 1:N_mcmc) {
    # 1. Actualizar h_t para cada t
    for (ti in 1:T) {
      res_h <- actualizar_h_t(t = ti, h = h, y = y, mu = mu, phi_h = phi_h,
                              sigma_h = sigma_h, proposal_sd = proposal_sd_h)
      h[ti] <- res_h$h_t_nuevo
      acept_h_total <- acept_h_total + res_h$aceptado
    }
    
    # 2. Actualizar conjuntamente (mu, phi_h) mediante Gibbs
    mu_phi_new <- actualizar_mu_phi_joint(h = h, sigma_h = sigma_h, b0 = b0, Sigma0 = Sigma0)
    # Si la propuesta genera phi_h fuera de (0,1), se podría rechazar o ajustar;
    # en este ejemplo, se asume que la contribución posterior es muy baja en esos casos.
    if (mu_phi_new[2] <= -2 || mu_phi_new[2] >= 2) {
      # Mantener valores actuales
      mu_new <- mu
      phi_new <- phi_h
    } else {
      mu_new <- mu_phi_new[1]
      phi_new <- mu_phi_new[2]
      acept_mu_phi <- acept_mu_phi + 1
    }
    mu <- mu_new
    phi_h <- phi_new
    
    # 3. Actualizar sigma_h
    sigma_h <- actualizar_sigma_h(h = h, mu = mu, phi_h = phi_h,
                                  prior_sigma_shape = prior_sigma_shape,
                                  prior_sigma_rate = prior_sigma_rate)
    
    # 4. Guardar muestras
    mu_save[iter] <- mu
    phi_h_save[iter] <- phi_h
    sigma_h_save[iter] <- sigma_h
    h_save[iter, ] <- h
    
    if (iter %% 1000 == 0) {
      cat("Iter:", iter, "mu:", round(mu,3), "phi_h:", round(phi_h,3),
          "sigma_h:", round(sigma_h,3), "\n")
    }
  }
  
  tasa_acept_mu_phi <- acept_mu_phi / N_mcmc
  tasa_acept_h <- acept_h_total / (N_mcmc * T)
  
  return(list(mu = mu_save, phi_h = phi_h_save, sigma_h = sigma_h_save, h = h_save,
              tasa_acept_mu_phi = tasa_acept_mu_phi, tasa_acept_h = tasa_acept_h))
}
