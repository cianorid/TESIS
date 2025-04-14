
filtro_particulas <- function(y_data, N_particulas, mu, phi_h, sigma_h) {
  T <- length(y_data)
  h_est <- numeric(T)
  
  h_particulas <- rnorm(N_particulas, mean = mu, sd = sigma_h)
  pesos <- rep(1/N_particulas, N_particulas)
  
  for (t in 1:T) {
    h_particulas <- mu + phi_h * (h_particulas - mu) + sigma_h * rnorm(N_particulas)
    
    log_pesos <- -0.5 * (y_data[t]^2 / exp(h_particulas) + h_particulas)
    pesos <- exp(log_pesos - max(log_pesos))  
    pesos <- pesos / sum(pesos)
    
    indices <- sample(1:N_particulas, size = N_particulas, replace = TRUE, prob = pesos)
    h_particulas <- h_particulas[indices]
    
    h_est[t] <- mean(h_particulas)
  }
  
  return(h_est)
}



