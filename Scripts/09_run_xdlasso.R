# ============================================================
# 09_run_xdlasso.R
# XD-Lasso — Phillips Curve: United States and Euro Area
# ============================================================

library(LasForecast)
library(dplyr)
library(tibble)

source("R/ivx_functions.R")
source("Scripts/01_build_ea_data.R")
source("Scripts/02_build_ea_data_quarterly.R")
source("Scripts/03_build_us_data.R")
source("Scripts/04_build_us_data_quarterly.R")

# ============================================================
# FUNCIÓN PRINCIPAL
# ============================================================

main_xdlasso <- function(data_md, y_name, x_name, include_vars,
                         sub_periods, vrtest = FALSE,
                         a_val = 0.5, se_type = "robust") {

  data_md <- data_md %>%
    select(all_of(c("date", y_name, x_name, include_vars)))

  results_list <- vector("list", length(sub_periods))

  for (i in seq_along(sub_periods)) {
    period_name  <- names(sub_periods)[i]
    period_dates <- sub_periods[[i]]

    data_sub <- data_md %>%
      filter(date >= period_dates[1] & date <= period_dates[2])

    T_sub <- nrow(data_sub)
    y <- data_sub[[y_name]]
    d <- data_sub[[x_name]]
    x <- as.matrix(data_sub %>% select(all_of(include_vars)))

    w_est <- cbind(d[-T_sub], x[-T_sub, ])
    y_est <- y[-1]
    d_est <- as.matrix(d[-T_sub])

    y_est <- y_est - mean(y_est, na.rm = TRUE)
    w_est <- apply(w_est, 2, function(col) col - mean(col, na.rm = TRUE))
    d_est <- apply(d_est, 2, function(col) col - mean(col, na.rm = TRUE))

    n_eff <- length(y_est)
    cat("  Periodo:", period_name, "| T =", n_eff, "| p =", ncol(w_est), "\n")

    deb_out <- debias_ivx(
      w            = w_est,
      y            = y_est,
      d_ind        = 1,
      train_method = "cv",
      k            = 10,
      a            = a_val,
      c_z          = 5,
      se_type      = se_type,
      zhang_zhang  = TRUE
    )

    vr_p <- if (vrtest) {
      vrtest::AutoBoot.test(deb_out$u_hat, nboot = 5000, wild = "Normal")$pval
    } else NA

    res_ivx <- ivx_inference(
      w   = as.matrix(d_est),
      y   = y_est,
      a   = a_val,
      c_z = 5
    )

    results_list[[i]] <- tibble(
      period      = period_name,
      n           = n_eff,
      ivx_est     = as.numeric(res_ivx$iv_est),
      ivx_se      = as.numeric(res_ivx$iv_se),
      dlasso_est  = deb_out$theta_hat_zz,
      dlasso_se   = deb_out$sigma_hat_zz,
      xdlasso_est = deb_out$theta_hat_ivx,
      xdlasso_se  = deb_out$sigma_hat_ivx,
      vrtest_pval = vr_p
    )
  }

  bind_rows(results_list)
}

# ============================================================
# FUNCIÓN DE IMPRESIÓN — tabla encajonada con asteriscos
# ============================================================

sig_stars <- function(t) {
  if (is.na(t)) return("   ")
  a <- abs(t)
  if (a > 2.576) return("***")
  if (a > 1.960) return(" **")
  if (a > 1.645) return("  *")
  return("   ")
}

print_xdlasso_table <- function(res, title, subtitle, show_vr = FALSE) {
  W_data <- if (show_vr) 109 else 103
  W <- W_data

  sep_e <- paste0("+", strrep("=", W), "+")
  sep_d <- paste0("+", strrep("-", W), "+")

  center <- function(s, w) {
    n <- nchar(s); p <- max(0, w - n)
    paste0(strrep(" ", floor(p/2)), s, strrep(" ", ceiling(p/2)))
  }
  rpad <- function(s, w) {
    n <- nchar(s)
    if (n < w) paste0(s, strrep(" ", w - n)) else substr(s, 1, w)
  }

  cat("\n")
  cat(sep_e, "\n", sep = "")
  cat(sprintf("| %s |\n", center(title,    W - 2)))
  cat(sprintf("| %s |\n", center(subtitle, W - 2)))
  cat(sep_e, "\n", sep = "")

  if (show_vr) {
    hdr <- sprintf("| %-15s | %5s | %8s %8s %4s | %8s %8s %4s | %8s %8s %4s | %5s |",
                   "Periodo", "T",
                   "IVX", "(SE)", "[t]",
                   "Dlasso", "(SE)", "[t]",
                   "XDlasso", "(SE)", "[t]", "VR_p")
  } else {
    hdr <- sprintf("| %-15s | %5s | %8s %8s %4s | %8s %8s %4s | %8s %8s %4s |",
                   "Periodo", "T",
                   "IVX", "(SE)", "[t]",
                   "Dlasso", "(SE)", "[t]",
                   "XDlasso", "(SE)", "[t]")
  }
  cat(hdr, "\n", sep = "")
  cat(sep_d, "\n", sep = "")

  for (i in seq_len(nrow(res))) {
    ivx_t <- res$ivx_est[i]     / res$ivx_se[i]
    dl_t  <- res$dlasso_est[i]  / res$dlasso_se[i]
    xdl_t <- res$xdlasso_est[i] / res$xdlasso_se[i]

    if (show_vr) {
      cat(sprintf("| %-15s | %5d | %8.3f %8s %4s | %8.3f %8s %4s | %8.3f %8s %4s | %5.3f |\n",
                  res$period[i], res$n[i],
                  res$ivx_est[i],     sprintf("(%6.3f)", res$ivx_se[i]),     sig_stars(ivx_t),
                  res$dlasso_est[i],  sprintf("(%6.3f)", res$dlasso_se[i]),  sig_stars(dl_t),
                  res$xdlasso_est[i], sprintf("(%6.3f)", res$xdlasso_se[i]), sig_stars(xdl_t),
                  res$vrtest_pval[i]))
    } else {
      cat(sprintf("| %-15s | %5d | %8.3f %8s %4s | %8.3f %8s %4s | %8.3f %8s %4s |\n",
                  res$period[i], res$n[i],
                  res$ivx_est[i],     sprintf("(%6.3f)", res$ivx_se[i]),     sig_stars(ivx_t),
                  res$dlasso_est[i],  sprintf("(%6.3f)", res$dlasso_se[i]),  sig_stars(dl_t),
                  res$xdlasso_est[i], sprintf("(%6.3f)", res$xdlasso_se[i]), sig_stars(xdl_t)))
    }
  }

  cat(sep_d, "\n", sep = "")
  nota1 <- "  *** p<0.01  ** p<0.05  * p<0.10  (|t-stat|)"
  nota2 <- "  (SE) entre parentesis  |  [t] = coef / SE"
  nota3 <- "  VR_p = p-valor Variance Ratio Test (Kim 2009)"
  cat(sprintf("| %s |\n", rpad(nota1, W - 2)))
  cat(sprintf("| %s |\n", rpad(nota2, W - 2)))
  if (show_vr) cat(sprintf("| %s |\n", rpad(nota3, W - 2)))
  cat(sep_e, "\n\n", sep = "")
}

# ============================================================
# VARIABLES LABORALES EA (para ejercicio Include / Exclude)
# ============================================================

labor_vars_ea <- c(
  "TEMP_EA", "EMP_EA", "SEMP_EA", "THOURS_EA",
  "EMPAG_EA", "EMPIN_EA", "EMPMN_EA", "EMPCON_EA",
  "EMPRT_EA", "EMPIT_EA", "EMPFC_EA", "EMPRE_EA",
  "EMPPR_EA", "EMPPA_EA", "EMPENT_EA", "RPRP_EA",
  "WS_EA", "ESC_EA",
  "ULCIN_EA", "ULCMQ_EA", "ULCMN_EA", "ULCCON_EA",
  "ULCRT_EA", "ULCFC_EA", "ULCRE_EA", "ULCPR_EA"
)

# ============================================================
# 1. USA — MENSUAL
# ============================================================

cat("\n====================================================\n")
cat(" 1. USA MENSUAL — XDlasso\n")
cat("====================================================\n\n")

sub_periods_us_m <- list(
  full_sample = c(as.Date("1960-01-01"), as.Date("2025-04-01")),
  pre_volcker = c(as.Date("1960-01-01"), as.Date("1979-07-01"))
)

cat("--- Panel transformado ---\n")
res_us_m_trans <- main_xdlasso(
  data_md = data_trans_us, y_name = "infl_cpi", x_name = "UNRATE_level",
  include_vars = vars_us, sub_periods = sub_periods_us_m,
  vrtest = FALSE, se_type = "robust", a_val = 0.5
)

cat("--- Panel sin transformar ---\n")
res_us_m_notrans <- main_xdlasso(
  data_md = data_notrans_us, y_name = "infl_cpi", x_name = "UNRATE_level",
  include_vars = vars_notrans_us, sub_periods = sub_periods_us_m,
  vrtest = FALSE, se_type = "robust", a_val = 0.5
)

print_xdlasso_table(
  res_us_m_trans,
  title    = "TABLA A1 — USA Mensual | Panel Transformado (TCODE)",
  subtitle = sprintf("Phillips: desempleo -> CPI | p = %d controles", length(vars_us)),
  show_vr  = FALSE
)

print_xdlasso_table(
  res_us_m_notrans,
  title    = "TABLA A2 — USA Mensual | Panel Sin Transformar (excl. I(2))",
  subtitle = sprintf("Phillips: desempleo -> CPI | p = %d controles", length(vars_notrans_us)),
  show_vr  = FALSE
)

# ============================================================
# 2. USA — TRIMESTRAL
# ============================================================

cat("\n====================================================\n")
cat(" 2. USA TRIMESTRAL — XDlasso\n")
cat("====================================================\n\n")

sub_periods_us_q <- list(
  full_sample = c(as.Date("1960-01-01"), as.Date("2025-07-01")),
  pre_volcker = c(as.Date("1960-01-01"), as.Date("1979-07-01"))
)

cat("--- Panel transformado ---\n")
res_us_q_trans <- main_xdlasso(
  data_md = data_trans_us_q, y_name = "infl_cpi", x_name = "UNRATE_level",
  include_vars = vars_us_q, sub_periods = sub_periods_us_q,
  vrtest = FALSE, se_type = "robust", a_val = 0.5
)

cat("--- Panel sin transformar ---\n")
res_us_q_notrans <- main_xdlasso(
  data_md = data_notrans_us_q, y_name = "infl_cpi", x_name = "UNRATE_level",
  include_vars = vars_notrans_us_q, sub_periods = sub_periods_us_q,
  vrtest = FALSE, se_type = "robust", a_val = 0.5
)

print_xdlasso_table(
  res_us_q_trans,
  title    = "TABLA B1 — USA Trimestral | Panel Transformado (TCODE)",
  subtitle = sprintf("Phillips: desempleo -> CPI | p = %d controles", length(vars_us_q)),
  show_vr  = FALSE
)

print_xdlasso_table(
  res_us_q_notrans,
  title    = "TABLA B2 — USA Trimestral | Panel Sin Transformar (excl. I(2))",
  subtitle = sprintf("Phillips: desempleo -> CPI | p = %d controles", length(vars_notrans_us_q)),
  show_vr  = FALSE
)

# ============================================================
# 3. EA — MENSUAL
# ============================================================

cat("\n====================================================\n")
cat(" 3. EA MENSUAL — XDlasso\n")
cat("====================================================\n\n")

sub_periods_ea_m <- list(
  full_sample = c(as.Date("2000-03-01"), as.Date("2025-02-01")),
  pre_2008    = c(as.Date("2000-03-01"), as.Date("2007-12-01")),
  post_2008   = c(as.Date("2008-01-01"), as.Date("2025-02-01"))
)

cat("--- Panel transformado | Include Labor ---\n")
res_ea_m_trans_incl <- main_xdlasso(
  data_md = data_trans_ea, y_name = "infl_hicp", x_name = "UNETOT_level",
  include_vars = vars_ea, sub_periods = sub_periods_ea_m,
  vrtest = FALSE, se_type = "robust", a_val = 0.5
)

cat("--- Panel sin transformar | Include Labor ---\n")
res_ea_m_notrans_incl <- main_xdlasso(
  data_md = data_notrans_ea, y_name = "infl_hicp", x_name = "UNETOT_level",
  include_vars = vars_notrans_ea, sub_periods = sub_periods_ea_m,
  vrtest = FALSE, se_type = "robust", a_val = 0.5
)

print_xdlasso_table(
  res_ea_m_trans_incl,
  title    = "TABLA C1a — EA Mensual | Transformado | Include Labor",
  subtitle = sprintf("Phillips: desempleo -> HICP | p = %d controles", length(vars_ea)),
  show_vr  = FALSE
)
print_xdlasso_table(
  res_ea_m_notrans_incl,
  title    = "TABLA C2a — EA Mensual | Sin Transformar | Include Labor",
  subtitle = sprintf("Phillips: desempleo -> HICP | p = %d controles", length(vars_notrans_ea)),
  show_vr  = FALSE
)

# ============================================================
# 4. EA — TRIMESTRAL
# ============================================================

cat("\n====================================================\n")
cat(" 4. EA TRIMESTRAL — XDlasso\n")
cat("====================================================\n\n")

vars_ea_q_nolabor         <- setdiff(vars_ea_q,         labor_vars_ea)
vars_notrans_ea_q_nolabor <- setdiff(vars_notrans_ea_q, labor_vars_ea)

fecha_inicio_q <- data_trans_ea_q$date[1]
fecha_fin_q    <- data_trans_ea_q$date[nrow(data_trans_ea_q)]

sub_periods_ea_q <- list(
  full_sample = c(fecha_inicio_q,       fecha_fin_q),
  pre_2008    = c(fecha_inicio_q,       as.Date("2007-10-01")),
  post_2008   = c(as.Date("2008-01-01"), fecha_fin_q)
)

cat("Muestra trimestral:", as.character(fecha_inicio_q), "a",
    as.character(fecha_fin_q), "| T =", nrow(data_trans_ea_q), "\n")
cat("p transformado   — Include:", length(vars_ea_q),
    "| Exclude:", length(vars_ea_q_nolabor), "\n")
cat("p sin transformar — Include:", length(vars_notrans_ea_q),
    "| Exclude:", length(vars_notrans_ea_q_nolabor), "\n\n")

cat("--- Panel transformado | Include Labor ---\n")
res_ea_q_trans_incl <- main_xdlasso(
  data_md = data_trans_ea_q, y_name = "infl_hicp", x_name = "UNETOT_level",
  include_vars = vars_ea_q, sub_periods = sub_periods_ea_q,
  vrtest = TRUE, se_type = "robust", a_val = 0.5
)

cat("--- Panel transformado | Exclude Labor ---\n")
res_ea_q_trans_excl <- main_xdlasso(
  data_md = data_trans_ea_q, y_name = "infl_hicp", x_name = "UNETOT_level",
  include_vars = vars_ea_q_nolabor, sub_periods = sub_periods_ea_q,
  vrtest = TRUE, se_type = "robust", a_val = 0.5
)

cat("--- Panel sin transformar | Include Labor ---\n")
res_ea_q_notrans_incl <- main_xdlasso(
  data_md = data_notrans_ea_q, y_name = "infl_hicp", x_name = "UNETOT_level",
  include_vars = vars_notrans_ea_q, sub_periods = sub_periods_ea_q,
  vrtest = TRUE, se_type = "robust", a_val = 0.5
)

cat("--- Panel sin transformar | Exclude Labor ---\n")
res_ea_q_notrans_excl <- main_xdlasso(
  data_md = data_notrans_ea_q, y_name = "infl_hicp", x_name = "UNETOT_level",
  include_vars = vars_notrans_ea_q_nolabor, sub_periods = sub_periods_ea_q,
  vrtest = TRUE, se_type = "robust", a_val = 0.5
)

print_xdlasso_table(
  res_ea_q_trans_incl,
  title    = "TABLA D1a — EA Trimestral | Transformado | Include Labor",
  subtitle = sprintf("Phillips: desempleo -> HICP | p = %d controles | VR Test incluido", length(vars_ea_q)),
  show_vr  = TRUE
)
print_xdlasso_table(
  res_ea_q_trans_excl,
  title    = "TABLA D1b — EA Trimestral | Transformado | Exclude Labor",
  subtitle = sprintf("Phillips: desempleo -> HICP | p = %d controles | VR Test incluido", length(vars_ea_q_nolabor)),
  show_vr  = TRUE
)
print_xdlasso_table(
  res_ea_q_notrans_incl,
  title    = "TABLA D2a — EA Trimestral | Sin Transformar | Include Labor",
  subtitle = sprintf("Phillips: desempleo -> HICP | p = %d controles | VR Test incluido", length(vars_notrans_ea_q)),
  show_vr  = TRUE
)
print_xdlasso_table(
  res_ea_q_notrans_excl,
  title    = "TABLA D2b — EA Trimestral | Sin Transformar | Exclude Labor",
  subtitle = sprintf("Phillips: desempleo -> HICP | p = %d controles | VR Test incluido", length(vars_notrans_ea_q_nolabor)),
  show_vr  = TRUE
)

cat("✅ 06_run_xdlasso.R completado\n")
