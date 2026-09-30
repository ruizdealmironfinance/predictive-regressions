# ============================================================
# 00_ivx_functions.R
# MÉTODO — no tocar para cambiar datos
#
# Funciones:
#   build_instrument()   → construye el IVX (convención D&R o Gao)
#   fit_ar_regressor()   → AR(p) sobre el regresor para RA-IVX
#   calc_ivx_raivx()     → estima IVX y RA-IVX dado (Y, X, Z, NU)
#   run_ivx_pipeline()   → orquesta todo desde un data.frame limpio
#
# CONTRATO DE INPUT — igual para todos los problemas:
#   df       → data.frame {date, y, x} longitud T_eff
#              y = variable dependiente ya transformada
#              x = regresor (puede ser recortado, solo para alineación)
#   x_full   → vector del regresor completo, longitud T_eff + 1
#              necesario para lags e instrumento IVX
#   convention → "DR" | "Gao"
#   label      → string para imprimir
#
# T_eff = nrow(df) siempre — sin casos especiales entre problemas
#
# NOTA SOBRE NU:
#   NU_dm        = nu_hat[t+1] — innovación futura, input para RA-IVX y gamma_hat
#   NU_contemp_dm = nu_hat[t]  — innovación contemporánea, solo para diagnósticos
#                                (cor_uv, delta_hat = parámetros de Stambaugh bias)
# ============================================================

library(sandwich)

# ------------------------------------------------------------
# 1. Instrumento IVX — construido sobre x_full completo
# ------------------------------------------------------------

build_instrument <- function(x_full, T_eff, convention = c("DR", "Gao")) {
  convention <- match.arg(convention)
  T_full <- length(x_full)   # T_eff + 1
  dx     <- diff(x_full)
  
  if (convention == "DR") {
    rho_z <- 1 - T_eff^(-0.95)
  } else {
    rho_z <- 1 - 5 / T_eff^0.5
  }
  
  z    <- numeric(T_full)
  z[1] <- 0
  for (t in 2:T_full) z[t] <- rho_z * z[t - 1] + dx[t - 1]
  
  list(z = z, rho_z = rho_z)
}

# ------------------------------------------------------------
# 2. AR(p) sobre x_full completo → residuos para RA-IVX
# ------------------------------------------------------------

fit_ar_regressor <- function(x_full, T_eff) {
  p_max  <- floor(4 * (T_eff / 100)^0.25)
  ar_fit <- ar(x_full, aic = TRUE, order.max = p_max,
               method = "ols", demean = FALSE)
  
  list(
    p       = ar_fit$order,
    nu_hat  = as.numeric(ar_fit$resid),
    rho_cap = round(sum(ar_fit$ar), 4),
    ar_fit  = ar_fit
  )
}

# ------------------------------------------------------------
# 3. Estimación IVX y RA-IVX
#
# CAMBIO: añadido argumento NU_contemp_dm
#   NU_dm         → nu_hat[t+1], demeaned — para gamma_hat y RA-IVX (sin cambio)
#   NU_contemp_dm → nu_hat[t],   demeaned — para cor_uv y delta_hat (Stambaugh)
# ------------------------------------------------------------

calc_ivx_raivx <- function(Y_dm, Xl_dm, NU_dm, Zl,
                           x_full, idx_eff, n_common, p,
                           NU_contemp_dm) {          # ← NUEVO argumento
  
  # -- IVX --
  beta_ivx         <- sum(Zl * Y_dm) / sum(Zl * Xl_dm)
  beta_ols_for_ivx <- sum(Xl_dm * Y_dm) / sum(Xl_dm^2)
  res_ols          <- Y_dm - beta_ols_for_ivx * Xl_dm
  denom_ivx        <- abs(sum(Zl * Xl_dm))
  sigma2_ivx       <- sum(res_ols^2) / n_common
  se_ivx           <- sqrt(sum(Zl^2 * res_ols^2)) / denom_ivx
  t_ivx            <- beta_ivx / se_ivx
  
  # -- RA-IVX --
  # gamma_hat usa NU_dm (t+1) — correcto, no cambia
  gamma_hat  <- sum(NU_dm * Y_dm) / sum(NU_dm^2)
  Y_tilde    <- Y_dm - gamma_hat * NU_dm
  beta_raivx <- sum(Zl * Y_tilde) / sum(Zl * Xl_dm)
  
  XZ           <- cbind(Xl_dm, NU_dm)
  ols_coefs_ra <- solve(t(XZ) %*% XZ) %*% (t(XZ) %*% Y_dm)
  eps_tilde    <- as.numeric(Y_dm - XZ %*% ols_coefs_ra)
  
  X_lags <- matrix(NA, nrow = n_common, ncol = p)
  for (j in 1:p) {
    vals        <- x_full[idx_eff - j + 1]
    X_lags[, j] <- vals - mean(vals)
  }
  
  H_xx_sum   <- t(X_lags) %*% X_lags
  H_zx_sum   <- t(Zl) %*% X_lags
  H_xxnu_sum <- t(X_lags) %*% (X_lags * NU_dm^2)
  
  Q_hat <- as.numeric(
    H_zx_sum %*% solve(H_xx_sum) %*% H_xxnu_sum %*%
      solve(H_xx_sum) %*% t(H_zx_sum)
  )
  
  denom_raivx      <- abs(sum(Zl * Xl_dm))
  EW_term          <- sum(Zl^2 * eps_tilde^2)
  gamma2Q          <- gamma_hat^2 * Q_hat
  num_var          <- EW_term + gamma2Q
  se_raivx         <- sqrt(num_var) / denom_raivx
  t_raivx          <- beta_raivx / se_raivx
  
  ratio_Q          <- gamma2Q / EW_term
  inflation_factor <- sqrt(num_var / EW_term)   # SE_raivx / SE_ivx
  
  # -- Stambaugh diagnostics — usan NU_contemp_dm (t), no NU_dm (t+1) --
  # CAMBIO: cor_uv, sigma_uv, delta_hat ahora usan NU_contemp_dm
  sigma2_u    <- sum(res_ols^2)       / n_common
  sigma2_v    <- sum(NU_contemp_dm^2) / n_common          # ← CAMBIADO
  sigma_uv    <- sum(res_ols * NU_contemp_dm) / n_common  # ← CAMBIADO
  delta_hat   <- sigma_uv / sigma2_v                       # = σ_uv / σ²_v
  cor_uv      <- cor(res_ols, NU_contemp_dm,               # ← CAMBIADO
                     use = "complete.obs")
  IF_stambaugh <- 1 + delta_hat * (sigma_uv / sigma2_u)
  
  list(
    beta_ivx         = beta_ivx,
    t_ivx            = t_ivx,
    beta_raivx       = beta_raivx,
    t_raivx          = t_raivx,
    gamma_hat        = gamma_hat,
    cor_uv           = cor_uv,
    EW_term          = EW_term,
    Q_hat            = Q_hat,
    gamma2Q          = gamma2Q,
    ratio_Q          = ratio_Q,
    inflation_factor = inflation_factor,
    IF_stambaugh     = IF_stambaugh,
    delta_hat        = delta_hat,
    denom_raivx      = denom_raivx
  )
}

# ------------------------------------------------------------
# 4. Pipeline completo
#
# ALINEACIÓN — igual en todos los problemas:
#
#   df$y     longitud T_eff  → variable dependiente ya transformada
#   df$x     longitud T_eff  → regresor alineado con y (solo OLS/NW)
#   x_full   longitud T_eff+1 → regresor completo para AR, IVX, lags
#
#   AR(p) sobre x_full → nu_hat longitud T_eff+1
#   Burn-in: primeros p valores de nu_hat son NA
#
#   idx_common = (p+1):T_eff
#   Y          = df$y[idx_common]
#   Xl         = df$x[idx_common]
#   Zl         = z_full[idx_common]
#   NU         = nu_hat[idx_common + 1]  ← innovación en t+1, para RA-IVX
#   NU_contemp = nu_hat[idx_common]      ← innovación en t,   para diagnósticos
#
#   n_common = T_eff - p  (tras complete.cases)
# ------------------------------------------------------------

run_ivx_pipeline <- function(df, x_full, convention = "DR", label = "") {
  
  stopifnot(all(c("date", "y", "x") %in% names(df)))
  stopifnot(length(x_full) == nrow(df) + 1)
  
  T_eff <- nrow(df)
  y_raw <- df$y
  dates <- df$date
  
  # -- AR sobre x_full completo (T_eff + 1 observaciones)
  ar_res <- fit_ar_regressor(x_full, T_eff = T_eff)
  p      <- ar_res$p
  nu_hat <- ar_res$nu_hat
  
  # -- Instrumento sobre x_full completo
  iv_res <- build_instrument(x_full, T_eff = T_eff, convention = convention)
  z_full <- iv_res$z
  
  # -- Alineación
  idx_common <- (p + 1):T_eff
  Y          <- y_raw[idx_common]
  Xl         <- df$x[idx_common]
  Zl         <- z_full[idx_common]
  NU         <- nu_hat[idx_common + 1]   # t+1 — para RA-IVX
  NU_contemp <- nu_hat[idx_common]       # t   — para diagnósticos Stambaugh  ← NUEVO
  
  # CAMBIO: valid incluye NU_contemp
  valid      <- complete.cases(Y, Xl, Zl, NU, NU_contemp)
  Y          <- Y[valid];  Xl <- Xl[valid]
  Zl         <- Zl[valid]; NU <- NU[valid]
  NU_contemp <- NU_contemp[valid]                                               # ← NUEVO
  idx_eff    <- idx_common[valid]
  n_common   <- length(Y)
  
  # -- Demeaning
  Y_dm          <- Y  - mean(Y)
  Xl_dm         <- Xl - mean(Xl)
  NU_dm         <- NU - mean(NU)
  NU_contemp_dm <- NU_contemp - mean(NU_contemp)                                # ← NUEVO
  
  # -- OLS + Newey-West
  lm_fit   <- lm(Y ~ Xl)
  beta_ols <- coef(lm_fit)["Xl"]
  t_ols    <- summary(lm_fit)$coefficients["Xl", "t value"]
  nw_vcov  <- sandwich::NeweyWest(lm_fit, lag = 1,
                                  prewhite = FALSE, adjust = FALSE)
  se_nw    <- sqrt(nw_vcov["Xl", "Xl"])
  t_nw     <- beta_ols / se_nw
  
  # -- IVX / RA-IVX — se pasa NU_contemp_dm como nuevo argumento
  res <- calc_ivx_raivx(Y_dm, Xl_dm, NU_dm, Zl,
                        x_full, idx_eff, n_common, p,
                        NU_contemp_dm)                                         # ← NUEVO
  
  # -- Tabla
  cat("\n=== RESULTADOS:", label, "| Convención:", convention, "===\n")
  cat("Muestra:", as.character(dates[idx_eff[1]]),
      "a", as.character(dates[idx_eff[n_common]]),
      "| n =", n_common, "\n")
  cat("AR(p): p =", p, "| rho_cap =", ar_res$rho_cap,
      "| rho_z =", round(iv_res$rho_z, 6), "\n\n")
  
  se_ols_explicit   <- summary(lm_fit)$coefficients["Xl", "Std. Error"]
  se_ivx_explicit   <- res$beta_ivx   / res$t_ivx
  se_raivx_explicit <- res$beta_raivx / res$t_raivx
  
  out <- data.frame(
    Estimador = c("OLS", "OLS-NW", "IVX", "RA-IVX"),
    Beta      = round(c(beta_ols, beta_ols,
                        res$beta_ivx, res$beta_raivx), 6),
    SE        = round(c(se_ols_explicit, se_nw,
                        se_ivx_explicit, se_raivx_explicit), 6),
    t_stat    = round(c(t_ols, t_nw,
                        res$t_ivx, res$t_raivx), 3),
    Sig_5pct  = abs(c(t_ols, t_nw,
                      res$t_ivx, res$t_raivx)) > 1.96
  )
  print(out, row.names = FALSE)
  cat("gamma_hat:", round(res$gamma_hat, 4),
      "| delta_hat:", round(res$delta_hat, 4),
      "| IF_stambaugh:", round(res$IF_stambaugh, 4),
      "| ratio_Q:", round(res$ratio_Q, 4),
      "| inflation_factor:", round(res$inflation_factor, 4),
      "| cor_uv:", round(res$cor_uv, 4), "\n")
  
  invisible(list(
    results  = out,
    res_full = res,
    
    diagnostics = list(
      label      = label,
      convention = convention,
      
      n          = n_common,
      date_start = dates[idx_eff[1]],
      date_end   = dates[idx_eff[n_common]],
      
      p          = p,
      rho_cap    = ar_res$rho_cap,
      rho_z      = iv_res$rho_z,
      
      gamma_hat        = res$gamma_hat,
      ratio_Q          = res$ratio_Q,
      inflation_factor = res$inflation_factor,
      cor_uv           = res$cor_uv,
      delta_hat        = res$delta_hat,
      IF_stambaugh     = res$IF_stambaugh
    ),
    
    ar       = ar_res,
    iv       = iv_res,
    idx_eff  = idx_eff,
    df_used  = df[idx_eff, ],
    
    n        = n_common,
    
    beta_ols = as.numeric(beta_ols),
    t_ols    = as.numeric(t_ols),
    t_nw     = as.numeric(t_nw)
  ))
}