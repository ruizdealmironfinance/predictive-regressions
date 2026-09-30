# ============================================================
# mc_ivx_simulation.R
# Monte Carlo: tamano empirico de OLS, IVX, RA-IVX
# ============================================================

source("R/ivx_functions.R")
library(MASS)

set.seed(2024)

fit_ar_regressor <- function(x_full, T_eff) {
  ar_fit <- ar(x_full, aic = FALSE, order.max = 1,
               method = "ols", demean = FALSE)
  list(
    p       = ar_fit$order,
    nu_hat  = as.numeric(ar_fit$resid),
    rho_cap = round(sum(ar_fit$ar), 4),
    ar_fit  = ar_fit
  )
}

make_df <- function(n, rho, rho_eu, beta = 0) {
  T_eff <- n
  Sigma  <- matrix(c(1, rho_eu, rho_eu, 1), 2, 2)
  shocks <- mvrnorm(T_eff + 2, mu = c(0, 0), Sigma = Sigma)
  eps    <- shocks[, 1]
  nu     <- shocks[, 2]
  
  x_full    <- numeric(T_eff + 1)
  x_full[1] <- 0
  for (t in seq_len(T_eff)) x_full[t + 1] <- rho * x_full[t] + nu[t + 1]
  
  y <- beta * x_full[1:T_eff] + eps[2:(T_eff + 1)]
  
  df <- data.frame(
    date = seq.Date(as.Date("1970-01-01"), by = "month", length.out = T_eff),
    y    = y,
    x    = x_full[1:T_eff]
  )
  list(df = df, x_full = x_full)
}

one_sim <- function(n, rho, rho_eu, beta = 0) {
  d <- make_df(n = n, rho = rho, rho_eu = rho_eu, beta = beta)
  capture.output(
    out <- run_ivx_pipeline(d$df, d$x_full, convention = "DR", label = "MC")
  )
  r <- out$results
  c(
    t_ols   = r$t_stat[r$Estimador == "OLS"],
    t_ivx   = r$t_stat[r$Estimador == "IVX"],
    t_raivx = r$t_stat[r$Estimador == "RA-IVX"]
  )
}

# ============================================================
# 3. MONTE CARLO
# ============================================================

N_sim    <- 5000
n_grid   <- c(50, 100, 200, 500)
rho_grid <- c(0.7, 0.9, 0.95, 0.99)
rho_eu   <- -0.90

cat(sprintf("Corriendo simulaciones... (%d replicas x %d celdas)\n",
            N_sim, length(n_grid) * length(rho_grid)))

results <- list()
k <- 1

for (n in n_grid) {
  for (rho in rho_grid) {
    t_mat <- t(replicate(N_sim, one_sim(n = n, rho = rho, rho_eu = rho_eu)))
    
    size_ols   <- mean(abs(t_mat[, "t_ols"])   > 1.96)
    size_ivx   <- mean(abs(t_mat[, "t_ivx"])   > 1.96)
    size_raivx <- mean(abs(t_mat[, "t_raivx"]) > 1.96)
    
    results[[k]] <- data.frame(
      n          = n,
      rho        = rho,
      size_OLS   = round(size_ols,   3),
      size_IVX   = round(size_ivx,   3),
      size_RAIVX = round(size_raivx, 3)
    )
    
    cat(sprintf("  n=%3d | rho=%.2f | OLS=%.3f | IVX=%.3f | RA-IVX=%.3f\n",
                n, rho, size_ols, size_ivx, size_raivx))
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

# ---- Exportar tabla de size ----
write.csv(
  res_df,
  "Results/Tables/Monte_Carlo/mc_size_results.csv",
  row.names = FALSE
)

# Tabla LaTeX para Overleaf
if (!requireNamespace("xtable", quietly = TRUE)) install.packages("xtable")
library(xtable)

size_tab <- res_df[order(res_df$n, res_df$rho), ]
colnames(size_tab) <- c("$n$", "$\\rho$", "OLS", "IVX", "RA-IVX")

print(
  xtable(size_tab,
         caption = "Empirical size at the 5\\% level ($\\beta=0$, $\\rho_{e\\nu}=-0.90$, $N_{sim}=2000$).",
         label   = "tab:mc_size",
         digits  = 3),
  include.rownames = FALSE,
  sanitize.colnames.function = identity,
  caption.placement = "top",
  file = "Results/Tables/Monte_Carlo/mc_size_table.tex"
)





# ============================================================
# 4. TABLA RESUMEN
# ============================================================

cat("\n")
cat(strrep("=", 65), "\n")
cat(sprintf("TAMANO EMPIRICO AL 5%% | beta=0, rho_eu=%.1f | N_sim=%d\n", rho_eu, N_sim))
cat(strrep("=", 65), "\n")
cat(sprintf("%-6s %-6s | %8s %8s %8s\n", "n", "rho", "OLS", "IVX", "RA-IVX"))
cat(strrep("-", 65), "\n")

for (i in seq_len(nrow(res_df))) {
  r <- res_df[i, ]
  cat(sprintf("%-6d %-6.2f | %8.3f %8.3f %8.3f\n",
              r$n, r$rho, r$size_OLS, r$size_IVX, r$size_RAIVX))
}

cat(strrep("=", 65), "\n")
cat(sprintf("Nominal: 0.050 | Banda 95%%: [0.036, 0.064] (N_sim=%d)\n", N_sim))
cat(strrep("=", 65), "\n\n")

# ============================================================
# 5. GRAFICO
# ============================================================

pdf("mc_size_plot.pdf", width = 11, height = 7)

par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 3, 0), bg = "white")

cols  <- c(OLS = "#e74c3c", IVX = "#2980b9", RAIVX = "#27ae60")
lty_v <- c(OLS = 2,         IVX = 1,          RAIVX = 1)
lwd_v <- c(OLS = 2,         IVX = 2.5,         RAIVX = 2.5)
pch_v <- c(OLS = 16,        IVX = 15,          RAIVX = 18)

for (n_val in n_grid) {
  sub   <- res_df[res_df$n == n_val, ]
  y_max <- max(sub$size_OLS, sub$size_IVX, sub$size_RAIVX, 0.20)
  
  plot(sub$rho, sub$size_OLS,
       type = "b", pch = pch_v["OLS"], col = cols["OLS"],
       lty = lty_v["OLS"], lwd = lwd_v["OLS"],
       ylim = c(0, y_max),
       xlab = expression(rho ~ "(predictor persistence)"),
       ylab = "Empirical size",
       main = bquote(n == .(n_val)),
       axes = FALSE)
  
  axis(1, at = rho_grid, labels = rho_grid)
  axis(2, las = 1)
  box()
  
  lines(sub$rho, sub$size_IVX,
        type = "b", pch = pch_v["IVX"], col = cols["IVX"],
        lty = lty_v["IVX"], lwd = lwd_v["IVX"])
  lines(sub$rho, sub$size_RAIVX,
        type = "b", pch = pch_v["RAIVX"], col = cols["RAIVX"],
        lty = lty_v["RAIVX"], lwd = lwd_v["RAIVX"])
  
  abline(h = 0.05,  col = "black",  lty = 2, lwd = 1)
  abline(h = 0.036, col = "grey60", lty = 3, lwd = 1)
  abline(h = 0.064, col = "grey60", lty = 3, lwd = 1)
  
  if (n_val == n_grid[1]) {
    legend("topleft",
           legend = c("OLS", "IVX", "RA-IVX", "Nominal 5%", "95% band"),
           col    = c(cols["OLS"], cols["IVX"], cols["RAIVX"], "black", "grey60"),
           lty    = c(lty_v["OLS"], lty_v["IVX"], lty_v["RAIVX"], 2, 3),
           lwd    = c(lwd_v["OLS"], lwd_v["IVX"], lwd_v["RAIVX"], 1, 1),
           pch    = c(pch_v["OLS"], pch_v["IVX"], pch_v["RAIVX"], NA, NA),
           bty    = "n", cex = 0.75)
  }
}

mtext(paste0("Empirical size at 5% | beta=0 | rho_eu=", rho_eu, " | N_sim=", N_sim),
      outer = TRUE, cex = 1.1, font = 2)

dev.off()
cat("Grafico guardado: mc_size_plot.pdf\n")
cat("mc_ivx_simulation.R completado\n")