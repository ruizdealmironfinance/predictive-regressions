# ============================================================
# 03_run_phillips_ivx.R
# SCRIPT MAESTRO — Curva de Phillips EA y USA
# IVX / RA-IVX univariados — Convención DR únicamente
# ============================================================

source("R/ivx_functions.R")

source("R/ivx_functions.R")

source("Scripts/01_build_ea_data.R")
source("Scripts/02_build_ea_data_quarterly.R")
source("Scripts/03_build_us_data.R")
source("Scripts/04_build_us_data_quarterly.R")

# ============================================================
# 1. ESTIMACIONES
# ============================================================

make_subperiod <- function(df, x_full, date_from, date_to) {
  idx <- which(df$date >= as.Date(date_from) & df$date <= as.Date(date_to))
  list(df = df[idx, ], x_full = x_full[c(idx, max(idx) + 1)])
}

# ── Área Euro — MENSUAL (full, como está en los datos) ───────
res_ea_dr        <- run_ivx_pipeline(df_ea, x_full_ea, convention = "DR", label = "EA mensual")

# ── Área Euro — MENSUAL subperiodos ──────────────────────────
sub_ea_pre08     <- make_subperiod(df_ea, x_full_ea, min(df_ea$date), "2007-12-01")
sub_ea_post08    <- make_subperiod(df_ea, x_full_ea, "2008-01-01",    max(df_ea$date))

res_ea_pre08_dr  <- run_ivx_pipeline(sub_ea_pre08$df,  sub_ea_pre08$x_full,  convention = "DR", label = "EA pre-2008")
res_ea_post08_dr <- run_ivx_pipeline(sub_ea_post08$df, sub_ea_post08$x_full, convention = "DR", label = "EA post-2008")

# ── Área Euro — TRIMESTRAL (full, como está en los datos) ────
res_ea_q_dr      <- run_ivx_pipeline(df_ea_q, x_full_ea_q, convention = "DR", label = "EA trimestral")

# ── USA — MENSUAL ────────────────────────────────────────────
sub_us_full_m    <- make_subperiod(df_us, x_full_us, "1959-01-01", "2025-04-01")
sub_pre79_m      <- make_subperiod(df_us, x_full_us, "1959-01-01", "1979-06-01")

res_us_dr        <- run_ivx_pipeline(sub_us_full_m$df, sub_us_full_m$x_full, convention = "DR", label = "USA mensual")
res_us_pre79_dr  <- run_ivx_pipeline(sub_pre79_m$df,   sub_pre79_m$x_full,   convention = "DR", label = "USA pre-Volcker")

# Great Moderation mensual: Stock & Watson (2002) — 1984M1-2007M12Q
sub_byp_m        <- make_subperiod(df_us, x_full_us, "1984-01-01", "2007-12-01")
res_us_byp_dr    <- run_ivx_pipeline(sub_byp_m$df, sub_byp_m$x_full, convention = "DR", label = "USA Great Moderation")

# ── USA — TRIMESTRAL ─────────────────────────────────────────
sub_us_full_q    <- make_subperiod(df_us_q, x_full_us_q, "1959-01-01", "2025-04-01")
sub_pre79_q      <- make_subperiod(df_us_q, x_full_us_q, "1959-01-01", "1979-07-01")

# Great Moderation trimestral
sub_byp_q        <- make_subperiod(df_us_q, x_full_us_q, "1984-01-01", "2007-10-01")

res_us_q_dr       <- run_ivx_pipeline(sub_us_full_q$df, sub_us_full_q$x_full, convention = "DR", label = "USA trimestral")
res_us_q_pre79_dr <- run_ivx_pipeline(sub_pre79_q$df,   sub_pre79_q$x_full,   convention = "DR", label = "USA pre-Volcker")
res_us_q_byp_dr   <- run_ivx_pipeline(sub_byp_q$df,     sub_byp_q$x_full,     convention = "DR", label = "USA Great Moderation")

# Volcker-Greenspan mensual: 1979M8 - 2006M1
sub_vg_m       <- make_subperiod(df_us, x_full_us, "1979-08-01", "2006-01-01")
res_us_vg_dr   <- run_ivx_pipeline(sub_vg_m$df, sub_vg_m$x_full, convention = "DR", label = "USA Volcker-Greenspan")

# Bernanke/Yellen/Powell mensual: 2006M2 - 2025M4
sub_bypow_m     <- make_subperiod(df_us, x_full_us, "2006-02-01", "2025-04-01")
res_us_bypow_dr <- run_ivx_pipeline(sub_bypow_m$df, sub_bypow_m$x_full, convention = "DR", label = "USA Bernanke/Yellen/Powell")

# Volcker-Greenspan trimestral: 1979Q4 - 2005Q4
# Volcker-Greenspan trimestral — mismas fechas exactas que mensual
sub_vg_q <- make_subperiod(df_us_q, x_full_us_q, "1979-07-01", "2006-01-01")
#sub_vg_q         <- make_subperiod(df_us_q, x_full_us_q, "1979-10-01", "2005-10-01")
res_us_q_vg_dr   <- run_ivx_pipeline(sub_vg_q$df, sub_vg_q$x_full, convention = "DR", label = "USA Volcker-Greenspan Q")

# Bernanke/Yellen/Powell trimestral: 2006Q1 - 2025Q2
sub_bypow_q       <- make_subperiod(df_us_q, x_full_us_q, "2006-01-01", "2025-04-01")
res_us_q_bypow_dr <- run_ivx_pipeline(sub_bypow_q$df, sub_bypow_q$x_full, convention = "DR", label = "USA Bernanke/Yellen/Powell Q")

# ============================================================
# 2. FUNCIONES AUXILIARES
# ============================================================

extract_row <- function(res, convencion) {
  r     <- res$results
  r_ols <- r[r$Estimador == "OLS",    ]
  r_nw  <- r[r$Estimador == "OLS-NW", ]
  r_ivx <- r[r$Estimador == "IVX",    ]
  r_ra  <- r[r$Estimador == "RA-IVX", ]
  
  data.frame(
    Convencion = convencion,
    
    beta_OLS = round(r_ols$Beta, 4),
    se_OLS   = round(r_ols$SE, 4),
    t_OLS    = round(r_ols$t_stat, 3),
    sig_OLS  = ifelse(abs(r_ols$t_stat) > 1.96, "*", " "),
    
    t_NW   = round(r_nw$t_stat, 3),
    se_NW  = round(r_nw$SE, 4),
    sig_NW = ifelse(abs(r_nw$t_stat) > 1.96, "*", " "),
    
    beta_IVX = round(r_ivx$Beta, 4),
    se_IVX   = round(r_ivx$SE, 4),
    t_IVX    = round(r_ivx$t_stat, 3),
    sig_IVX  = ifelse(abs(r_ivx$t_stat) > 1.96, "*", " "),
    
    beta_RAIVX = round(r_ra$Beta, 4),
    se_RAIVX   = round(r_ra$SE, 4),
    t_RAIVX    = round(r_ra$t_stat, 3),
    sig_RAIVX  = ifelse(abs(r_ra$t_stat) > 1.96, "*", " "),
    
    stringsAsFactors = FALSE
  )
}


extract_diag <- function(res) {
  d <- res$diagnostics
  
  # AR(1) simple del predictor (desempleo)
  x     <- res$df_used$x
  n_x   <- length(x)
  ar1_u <- round(coef(lm(x[-1] ~ x[-n_x]))[2], 4)
  
  data.frame(
    Muestra          = d$label,
    Convencion       = d$convention,
    n                = d$n,
    p                = d$p,
    ar1_u            = ar1_u,        # <-- sustituye rho_cap
    rho_z            = round(d$rho_z, 6),
    cor_uv           = round(d$cor_uv, 4),
    gamma_hat        = round(d$gamma_hat, 4),
    ratio_Q          = signif(d$ratio_Q, 3),
    delta_hat        = round(d$delta_hat, 4),
    IF_stambaugh     = round(d$IF_stambaugh, 4),
    inflation_factor = round(d$inflation_factor, 4),
    stringsAsFactors = FALSE
  )
}

print_ivx_table_old <- function(rows, title, subtitle) {
  W <- 110
  sep_e <- paste0("+", strrep("=", W), "+")
  sep_d <- paste0("+", strrep("-", W), "+")

  center <- function(s, w) {
    n <- nchar(s)
    p <- max(0, w - n)
    paste0(strrep(" ", floor(p / 2)), s, strrep(" ", ceiling(p / 2)))
  }

  rpad <- function(s, w) {
    n <- nchar(s)
    if (n < w) paste0(s, strrep(" ", w - n)) else substr(s, 1, w)
  }

  cat("\n")
  cat(sep_e, "\n", sep = "")
  cat(sprintf("| %s |\n", center(title, W - 2)))
  cat(sprintf("| %s |\n", center(subtitle, W - 2)))
  cat(sep_e, "\n", sep = "")

  cat(sprintf(
    "| %s | %s | %s |\n",
    center("", 10),
    center("--- OLS (sesgado) ---", 27),
    center("--- Estimadores Instrumentales ---", 68)
  ))

  cat(sprintf(
    "| %-10s | %8s | %7s | %1s | %7s | %1s || %8s | %7s | %1s | %8s | %7s | %1s |\n",
    "Convencion",
    "beta_OLS", "t_OLS", "*", "t_NW", "*",
    "beta_IVX", "t_IVX", "*", "beta_RAIVX", "t_RAIVX", "*"
  ))

  cat(sep_d, "\n", sep = "")

  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]

    cat(sprintf(
      "| %-10s | %8.4f | %7.3f | %1s | %7.3f | %1s || %8.4f | %7.3f | %1s | %8.4f | %7.3f | %1s |\n",
      r$Convencion,
      r$beta_OLS, r$t_OLS, r$sig_OLS,
      r$t_NW, r$sig_NW,
      r$beta_IVX, r$t_IVX, r$sig_IVX,
      r$beta_RAIVX, r$t_RAIVX, r$sig_RAIVX
    ))
  }

  cat(sep_d, "\n", sep = "")
  cat(sprintf("| %s |\n", rpad("  * significativo al 5% (|t| > 1.96)", W - 2)))
  cat(sprintf("| %s |\n", rpad("  OLS y OLS-NW comparten beta (mismo estimador puntual, distinto SE)", W - 2)))
  cat(sprintf("| %s |\n", rpad("  t_NW = t con errores Newey-West  ||  separa bloque OLS de IV", W - 2)))
  cat(sep_e, "\n\n", sep = "")
}
print_ivx_table_old1 <- function(rows, title, subtitle = NULL) {
  W <- 98
  sep_e <- paste0("+", strrep("=", W), "+")
  sep_d <- paste0("+", strrep("-", W), "+")
  
  center <- function(s, w) {
    s <- as.character(s)
    p <- max(0, w - nchar(s))
    paste0(strrep(" ", floor(p / 2)), s, strrep(" ", ceiling(p / 2)))
  }
  
  sig <- function(t) {
    ifelse(abs(t) > 2.576, "**",
           ifelse(abs(t) > 1.96, "*", ""))
  }
  
  fmt_num <- function(x, digits = 3) {
    sprintf(paste0("%.", digits, "f"), x)
  }
  
  fmt_t <- function(x) {
    paste0(fmt_num(x, 3), sig(x))
  }
  
  cat("\n")
  cat(sep_e, "\n", sep = "")
  cat(sprintf("| %s |\n", center(title, W - 2)))
  if (!is.null(subtitle)) {
    cat(sprintf("| %s |\n", center(subtitle, W - 2)))
  }
  cat(sep_e, "\n", sep = "")
  
  cat(sprintf(
    "| %-26s | %8s | %8s | %8s | %8s | %10s | %10s |\n",
    "Sample", "b_OLS", "t_OLS", "b_IVX", "t_IVX", "b_RA-IVX", "t_RA-IVX"
  ))
  
  cat(sep_d, "\n", sep = "")
  
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    
    cat(sprintf(
      "| %-26s | %8.4f | %8s | %8.4f | %8s | %10.4f | %10s |\n",
      r$Sample,
      r$beta_OLS,
      fmt_t(r$t_OLS),
      r$beta_IVX,
      fmt_t(r$t_IVX),
      r$beta_RAIVX,
      fmt_t(r$t_RAIVX)
    ))
  }
  
  cat(sep_d, "\n", sep = "")
  cat(sprintf("| %-96s |\n", "* p < 0.05, ** p < 0.01. OLS-NW omitted from this compact display."))
  cat(sep_e, "\n\n", sep = "")
}

make_row <- function(sample_name, res) {
  out <- extract_row(res, "DR")
  out$Sample <- sample_name
  
  out <- out[, c(
    "Sample",
    "beta_OLS", "se_OLS", "t_OLS",
    "beta_IVX", "se_IVX", "t_IVX",
    "beta_RAIVX", "se_RAIVX", "t_RAIVX"
  )]
  
  out
}

fmt_sig <- function(t) {
  ifelse(abs(t) > 2.576, "**",
         ifelse(abs(t) > 1.96, "*", ""))
}

fmt_num <- function(x, digits = 4) {
  sprintf(paste0("%.", digits, "f"), as.numeric(x))
}

fmt_t <- function(x) {
  paste0(sprintf("%.3f", as.numeric(x)), fmt_sig(as.numeric(x)))
}

fmt_beta_se <- function(beta, se) {
  sprintf("%.4f (%.4f)", as.numeric(beta), as.numeric(se))
}

print_panel <- function(title, rows) {
  
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  
  cat("\n")
  cat(strrep("=", 130), "\n")
  cat(title, "\n")
  cat(strrep("=", 130), "\n")
  
  cat(sprintf(
    "%-28s | %48s || %29s\n",
    "",
    "Coefficients (SE)",
    "t-statistics"
  ))
  
  cat(sprintf(
    "%-28s | %16s %16s %16s || %9s %9s %9s\n",
    "Sample",
    "OLS", "IVX", "RA-IVX",
    "OLS", "IVX", "RA-IVX"
  ))
  
  cat(strrep("-", 130), "\n")
  
  for(i in seq_len(nrow(rows))) {
    
    r <- rows[i, ]
    
    cat(sprintf(
      "%-28s | %16s %16s %16s || %9s %9s %9s\n",
      as.character(r$Sample),
      
      fmt_beta_se(r$beta_OLS, r$se_OLS),
      fmt_beta_se(r$beta_IVX, r$se_IVX),
      fmt_beta_se(r$beta_RAIVX, r$se_RAIVX),
      
      fmt_t(r$t_OLS),
      fmt_t(r$t_IVX),
      fmt_t(r$t_RAIVX)
    ))
  }
  
  cat(strrep("=", 130), "\n")
}


# ============================================================
# 3. CONSTRUIR FILAS DE LAS TABLAS
# ============================================================

# Panel A: Euro Area — Monthly
ea_monthly_rows <- rbind(
  make_row("Full sample",  res_ea_dr),
  make_row("Pre-2008",     res_ea_pre08_dr),
  make_row("Post-2008",    res_ea_post08_dr)
)

# Panel B: Euro Area — Quarterly
ea_quarterly_rows <- rbind(
  make_row("Full sample",  res_ea_q_dr)
)

# Panel C: United States — Monthly
us_monthly_rows <- rbind(
  make_row("Full sample",          res_us_dr),
  make_row("Pre-Volcker",          res_us_pre79_dr),
  make_row("Great Moderation",     res_us_byp_dr),
  make_row("Volcker-Greenspan",    res_us_vg_dr),
  make_row("Bernanke/Yellen/Powell", res_us_bypow_dr)
)

# Panel D: United States — Quarterly
us_quarterly_rows <- rbind(
  make_row("Full sample",          res_us_q_dr),
  make_row("Pre-Volcker",          res_us_q_pre79_dr),
  make_row("Great Moderation",     res_us_q_byp_dr),
  make_row("Volcker-Greenspan",    res_us_q_vg_dr),
  make_row("Bernanke/Yellen/Powell", res_us_q_bypow_dr)
)






print_panel("Panel A: Euro Area — Monthly", ea_monthly_rows)
print_panel("Panel B: Euro Area — Quarterly", ea_quarterly_rows)
print_panel("Panel C: United States — Monthly", us_monthly_rows)
print_panel("Panel D: United States — Quarterly", us_quarterly_rows)


# ============================================================
# EXPORTAR RESULTADOS PRINCIPALES
# ============================================================

ea_monthly_rows$Panel   <- "Euro Area - Monthly"
ea_quarterly_rows$Panel <- "Euro Area - Quarterly"
us_monthly_rows$Panel   <- "United States - Monthly"
us_quarterly_rows$Panel <- "United States - Quarterly"

main_results <- rbind(
  ea_monthly_rows,
  ea_quarterly_rows,
  us_monthly_rows,
  us_quarterly_rows
)

# Colocar Panel y Sample al principio
main_results <- main_results[, c(
  "Panel", "Sample",
  "beta_OLS", "se_OLS", "t_OLS",
  "beta_IVX", "se_IVX", "t_IVX",
  "beta_RAIVX", "se_RAIVX", "t_RAIVX"
)]

write.csv(
  main_results,
  "Results/Tables/Main/phillips_curve_results.csv",
  row.names = FALSE
)




# ============================================================
# 4. TABLA FINAL DE DIAGNÓSTICOS
# ============================================================

diag_table <- rbind(
  extract_diag(res_ea_dr),
  extract_diag(res_ea_pre08_dr),
  extract_diag(res_ea_post08_dr),
  extract_diag(res_ea_q_dr),
  extract_diag(res_us_dr),
  extract_diag(res_us_pre79_dr),
  extract_diag(res_us_byp_dr),
  extract_diag(res_us_q_dr),
  extract_diag(res_us_q_pre79_dr),
  extract_diag(res_us_q_byp_dr),
  extract_diag(res_us_vg_dr),
  extract_diag(res_us_bypow_dr),
  extract_diag(res_us_q_vg_dr),
  extract_diag(res_us_q_bypow_dr)
)

cat("\n")
cat("============================================================\n")
cat("DIAGNÓSTICOS — Persistencia y Endogeneidad\n")
cat("============================================================\n")
print(diag_table, row.names = FALSE)
# Exportar diagnósticos
dir.create(
  "Results/Tables/Diagnostics",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  diag_table,
  "Results/Tables/Diagnostics/ivx_diagnostics.csv",
  row.names = FALSE
)
cat("============================================================\n")
cat("rho_cap = suma de coeficientes AR(p) del predictor\n")
cat("rho_z   = parámetro de persistencia del instrumento IVX\n")
cat("cor_uv  = correlación entre residuo predictivo y shock del predictor\n")
cat("============================================================\n\n")

cat("Convencion DR = instrumento Demetrescu & Rodrigues (2022)\n")
cat("beta_IVX / beta_RAIVX = coeficiente estimado\n")
cat("t_IVX    / t_RAIVX    = estadistico t robusto\n\n")
cat("✅ 03_run_phillips_ivx.R completado\n")
