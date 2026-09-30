# ============================================================
# mc_power_rho99.R
# Potencia simetrica — solo rho = 0.99
# 6 paneles (uno por n), dentro de cada uno: OLS vs IVX vs RA-IVX
# ============================================================

source("R/ivx_functions.R")
library(MASS)

set.seed(2027)

fit_ar_regressor <- function(x_full, T_eff) {
  ar_fit <- ar(x_full, aic = FALSE, order.max = 1,
               method = "ols", demean = FALSE)
  list(p       = ar_fit$order,
       nu_hat  = as.numeric(ar_fit$resid),
       rho_cap = round(sum(ar_fit$ar), 4),
       ar_fit  = ar_fit)
}

make_df <- function(n, rho, rho_eu, beta = 0) {
  T_eff  <- n
  Sigma  <- matrix(c(1, rho_eu, rho_eu, 1), 2, 2)
  shocks <- mvrnorm(T_eff + 2, mu = c(0, 0), Sigma = Sigma)
  eps <- shocks[, 1]; nu <- shocks[, 2]
  x_full <- numeric(T_eff + 1); x_full[1] <- 0
  for (t in seq_len(T_eff)) x_full[t+1] <- rho * x_full[t] + nu[t+1]
  y <- beta * x_full[1:T_eff] + eps[2:(T_eff+1)]
  df <- data.frame(
    date = seq.Date(as.Date("1970-01-01"), by = "month", length.out = T_eff),
    y = y, x = x_full[1:T_eff])
  list(df = df, x_full = x_full)
}

one_sim <- function(n, rho, rho_eu, beta) {
  d <- make_df(n=n, rho=rho, rho_eu=rho_eu, beta=beta)
  capture.output(out <- run_ivx_pipeline(d$df, d$x_full, convention="DR", label="MC"))
  r <- out$results
  c(t_ols   = r$t_stat[r$Estimador == "OLS"],
    t_ivx   = r$t_stat[r$Estimador == "IVX"],
    t_raivx = r$t_stat[r$Estimador == "RA-IVX"])
}

N_sim     <- 2000
rho       <- 0.99
rho_eu    <- -0.9
n_grid    <- c(50, 100, 150, 200, 300, 500)
beta_grid <- seq(-0.15, 0.15, by = 0.025)

cat(sprintf("Corriendo... %d celdas x %d replicas\n",
            length(n_grid)*length(beta_grid), N_sim))

results <- list(); k <- 1
for (n in n_grid) {
  for (beta in beta_grid) {
    t_mat <- t(replicate(N_sim, one_sim(n=n, rho=rho, rho_eu=rho_eu, beta=beta)))
    results[[k]] <- data.frame(
      n         = n, beta = beta,
      pow_OLS   = mean(abs(t_mat[,"t_ols"])   > 1.96),
      pow_IVX   = mean(abs(t_mat[,"t_ivx"])   > 1.96),
      pow_RAIVX = mean(abs(t_mat[,"t_raivx"]) > 1.96))
    cat(sprintf("  n=%3d | beta=%+.3f | OLS=%.3f | IVX=%.3f | RA-IVX=%.3f\n",
                n, beta, results[[k]]$pow_OLS,
                results[[k]]$pow_IVX, results[[k]]$pow_RAIVX))
    k <- k + 1
  }
}
res_df <- do.call(rbind, results)




# Crear carpetas de resultados si no existen
dir.create(
  "Results/Tables/Monte_Carlo",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Figures/Monte_Carlo",
  recursive = TRUE,
  showWarnings = FALSE
)
# ---- Exportar tabla de power ----
write.csv(
  res_df,
  "Results/Tables/Monte_Carlo/mc_power_results.csv",
  row.names = FALSE
)


# ============================================================
# GRAFICO — 6 paneles (2 filas x 3 cols), uno por n
#           dentro de cada panel: OLS vs IVX vs RA-IVX
# ============================================================

est_cols <- c(OLS="#e74c3c", IVX="#2980b9", RAIVX="#27ae60")
est_lty  <- c(OLS=2,         IVX=1,          RAIVX=1)
est_lwd  <- c(OLS=1.8,       IVX=2,           RAIVX=2.5)
est_pch  <- c(OLS=16,        IVX=15,          RAIVX=18)

pdf(
  "Results/Figures/Monte_Carlo/mc_power_rho99.pdf",
  width = 13,
  height = 8
)
par(mfrow = c(2, 3), mar = c(4.5, 4, 3.5, 1),
    oma = c(0, 0, 3.5, 0), bg = "white")

for (n_val in n_grid) {

  sub <- res_df[res_df$n == n_val, ]
  sub <- sub[order(sub$beta), ]

  plot(NA, NA,
       xlim = c(-0.15, 0.15), ylim = c(0, 1),
       xlab = expression(beta ~ "(signal strength)"),
       ylab = "Power",
       main = bquote(n == .(n_val)),
       axes = FALSE)

  axis(1, at = beta_grid, labels = sprintf("%+.3f", beta_grid),
       cex.axis = 0.62, las = 2)
  axis(2, las = 1)
  box()

  abline(v = 0,    col = "grey70", lty = 1, lwd = 1)
  abline(h = 0.05, col = "grey50", lty = 2, lwd = 1)

  for (est in c("OLS", "IVX", "RAIVX")) {
    col_name <- paste0("pow_", est)
    lines(sub$beta, sub[[col_name]],
          col = est_cols[est], lty = est_lty[est], lwd = est_lwd[est])
    points(sub$beta, sub[[col_name]],
           col = est_cols[est], pch = est_pch[est], cex = 0.75)
  }

  if (n_val == n_grid[1]) {
    legend("topleft",
           legend = c("OLS", "IVX", "RA-IVX"),
           col    = est_cols,
           lty    = est_lty,
           lwd    = est_lwd,
           pch    = est_pch,
           bty    = "n", cex = 0.82)
  }
}

mtext(paste0("Power: OLS vs IVX vs RA-IVX | rho=0.99 | rho_eu=-0.9 | N_sim=", N_sim),
      outer = TRUE, cex = 1.05, font = 2)

dev.off()
cat("Guardado: mc_power_rho99.pdf\n")
cat("mc_power_rho99.R completado\n")

# ============================================================
# GRAFICO 2: Distribuciones del t-stat
# n=100, rho=0.99, tres paneles: beta=0, beta=+0.075, beta=-0.075
# ============================================================

cat("Generando distribuciones del t-stat...\n")

N_dist  <- 5000
n_dist  <- 100
rho_dist <- 0.99
betas_dist <- c(0, 0.075, -0.075)

# Guardar t-stats crudos (no solo rechazo)
one_sim_raw <- function(n, rho, rho_eu, beta) {
  d <- make_df(n=n, rho=rho, rho_eu=rho_eu, beta=beta)
  capture.output(
    out <- run_ivx_pipeline(d$df, d$x_full, convention="DR", label="MC")
  )
  r <- out$results
  c(t_ols   = r$t_stat[r$Estimador == "OLS"],
    t_ivx   = r$t_stat[r$Estimador == "IVX"],
    t_raivx = r$t_stat[r$Estimador == "RA-IVX"])
}

tdist_list <- list()
for (beta in betas_dist) {
  cat(sprintf("  Simulando beta=%+.3f...\n", beta))
  mat <- t(replicate(N_dist, one_sim_raw(n=n_dist, rho=rho_dist,
                                          rho_eu=-0.9, beta=beta)))
  tdist_list[[as.character(beta)]] <- mat
}

# --- Grafico ---
pdf(
  "Results/Figures/Monte_Carlo/mc_tdist_rho99.pdf",
  width = 13,
  height = 5
)
par(mfrow = c(1, 3), mar = c(4.5, 4, 3.5, 1),
    oma = c(0, 0, 3.5, 0), bg = "white")

est_cols <- c(OLS="#e74c3c", IVX="#2980b9", RAIVX="#27ae60")
est_lty  <- c(OLS=2, IVX=1, RAIVX=1)
est_lwd  <- c(OLS=1.8, IVX=2, RAIVX=2.5)

panel_titles <- c("0"    = expression(beta == 0 ~ "(H"[0]*")"),
                  "0.075"  = expression(beta == +0.075 ~ "(H"[1]*")"),
                  "-0.075" = expression(beta == -0.075 ~ "(H"[1]*")"))

for (beta in betas_dist) {
  mat  <- tdist_list[[as.character(beta)]]

  # Rango comun
  all_t <- c(mat[,"t_ols"], mat[,"t_ivx"], mat[,"t_raivx"])
  xlim  <- quantile(all_t, c(0.005, 0.995))
  xlim  <- c(max(xlim[1], -8), min(xlim[2], 8))

  # Densidades
  d_ols   <- density(mat[,"t_ols"],   bw="SJ", from=xlim[1], to=xlim[2])
  d_ivx   <- density(mat[,"t_ivx"],   bw="SJ", from=xlim[1], to=xlim[2])
  d_raivx <- density(mat[,"t_raivx"], bw="SJ", from=xlim[1], to=xlim[2])

  ylim <- c(0, max(d_ols$y, d_ivx$y, d_raivx$y) * 1.15)

  plot(NA, NA, xlim=xlim, ylim=ylim,
       xlab="t-statistic", ylab="Density",
       main=panel_titles[as.character(beta)],
       axes=TRUE)

  # Curva normal estandar de referencia
  x_ref <- seq(xlim[1], xlim[2], length=300)
  lines(x_ref, dnorm(x_ref), col="grey70", lty=3, lwd=1.5)

  # Zonas de rechazo
  abline(v =  1.96, col="grey50", lty=2, lwd=1)
  abline(v = -1.96, col="grey50", lty=2, lwd=1)

  # Densidades
  lines(d_ols$x,   d_ols$y,   col=est_cols["OLS"],
        lty=est_lty["OLS"],   lwd=est_lwd["OLS"])
  lines(d_ivx$x,   d_ivx$y,   col=est_cols["IVX"],
        lty=est_lty["IVX"],   lwd=est_lwd["IVX"])
  lines(d_raivx$x, d_raivx$y, col=est_cols["RAIVX"],
        lty=est_lty["RAIVX"], lwd=est_lwd["RAIVX"])

  # Medias verticales
  abline(v=mean(mat[,"t_ols"]),   col=est_cols["OLS"],
         lty=1, lwd=1.2)
  abline(v=mean(mat[,"t_ivx"]),   col=est_cols["IVX"],
         lty=1, lwd=1.2)
  abline(v=mean(mat[,"t_raivx"]), col=est_cols["RAIVX"],
         lty=1, lwd=1.2)

  if (beta == 0) {
    legend("topright",
           legend = c("OLS", "IVX", "RA-IVX", "N(0,1)"),
           col    = c(est_cols["OLS"], est_cols["IVX"],
                      est_cols["RAIVX"], "grey70"),
           lty    = c(est_lty["OLS"], est_lty["IVX"], 1, 3),
           lwd    = c(est_lwd, 1.5),
           bty    = "n", cex = 0.82)
  }

  # Anotar medias
  mtext(sprintf("mean: OLS=%.2f  IVX=%.2f  RA-IVX=%.2f",
                mean(mat[,"t_ols"]),
                mean(mat[,"t_ivx"]),
                mean(mat[,"t_raivx"])),
        side=3, line=0.2, cex=0.62, col="grey30")
}

mtext(paste0("t-stat distributions | n=", n_dist,
             " | rho=", rho_dist,
             " | rho_eu=-0.9 | N_sim=", N_dist),
      outer=TRUE, cex=1.0, font=2)

dev.off()
cat("Guardado: mc_tdist_rho99.pdf\n")
cat("Script completado\n")
