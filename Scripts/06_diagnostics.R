
# ============================================================
# 06_diagnostics.R
# Persistence diagnostics — Euro Area and United States
# ============================================================

source("R/ivx_functions.R")

source("Scripts/01_build_ea_data.R")
source("Scripts/02_build_ea_data_quarterly.R")
source("Scripts/03_build_us_data.R")
source("Scripts/04_build_us_data_quarterly.R")

library(urca)

# ------------------------------------------------------------
# Función
# ------------------------------------------------------------

persistence_row <- function(x) {
  x       <- as.numeric(na.omit(x))
  n       <- length(x)
  lag_adf <- floor(n^(1/3))
  
  # AR(1)
  ar1 <- round(coef(lm(x[-1] ~ x[-n]))[2], 3)
  
  # ADF con drift
  adf_res <- ur.df(x, type = "drift", lags = lag_adf)
  stat    <- as.numeric(adf_res@teststat[1])
  cv      <- adf_res@cval[1, ]
  
  # p-valor por tramos
  pval_s <- ifelse(stat < cv["1pct"],  "<0.01",
                   ifelse(stat < cv["5pct"],  "<0.05",
                          ifelse(stat < cv["10pct"], "<0.10",
                                 as.character(round(
                                   0.15 + (stat - cv["10pct"]) / abs(cv["10pct"]) * 0.15,
                                   3)))))
  
  list(ar1 = ar1, stat = round(stat, 3), pval = pval_s)
}

sub_x <- function(df, d1, d2) df$x[df$date >= as.Date(d1) & df$date <= as.Date(d2)]
sub_y <- function(df, d1, d2) df$y[df$date >= as.Date(d1) & df$date <= as.Date(d2)]

# ------------------------------------------------------------
# Muestras
# ------------------------------------------------------------

muestras <- list(
  # ── Euro Area ──────────────────────────────────────────────────────────────
  list(label="EA mensual (full)",   periodo="full sample", freq="M",
       y=df_ea$y,   x=df_ea$x),
  list(label="EA pre-2008",         periodo="hasta Dic 2007", freq="M",
       y=sub_y(df_ea,"1900-01-01","2007-12-01"),
       x=sub_x(df_ea,"1900-01-01","2007-12-01")),
  list(label="EA post-2008",        periodo="Ene 2008 – end", freq="M",
       y=sub_y(df_ea,"2008-01-01","2100-01-01"),
       x=sub_x(df_ea,"2008-01-01","2100-01-01")),
  list(label="EA trimestral (full)",periodo="full sample", freq="Q",
       y=df_ea_q$y, x=df_ea_q$x),
  
  # ── USA mensual ────────────────────────────────────────────────────────────
  list(label="USA mensual (full)",  periodo="Ene 1959 – Abr 2025", freq="M",
       y=df_us$y,   x=df_us$x),
  list(label="USA pre-Volcker",     periodo="Ene 1959 – Jun 1979", freq="M",
       y=sub_y(df_us,"1959-01-01","1979-06-01"),
       x=sub_x(df_us,"1959-01-01","1979-06-01")),
  list(label="USA Gran Moderación", periodo="Ene 1984 – Dic 2007", freq="M",
       y=sub_y(df_us,"1984-01-01","2007-12-01"),
       x=sub_x(df_us,"1984-01-01","2007-12-01")),
  list(label="USA Volcker-Greenspan",  periodo="Ago 1979 – Ene 2006", freq="M",
       y=sub_y(df_us,"1979-08-01","2006-01-01"),
       x=sub_x(df_us,"1979-08-01","2006-01-01")),
  list(label="USA Bernanke/Yellen/Powell", periodo="Feb 2006 – Abr 2025", freq="M",
       y=sub_y(df_us,"2006-02-01","2025-04-01"),
       x=sub_x(df_us,"2006-02-01","2025-04-01")),
  
  # ── USA trimestral ─────────────────────────────────────────────────────────
  list(label="USA trimestral (full)",      periodo="Ene 1959 – Abr 2025", freq="Q",
       y=df_us_q$y, x=df_us_q$x),
  list(label="USA pre-Volcker (trim.)",    periodo="Ene 1959 – Jul 1979", freq="Q",
       y=sub_y(df_us_q,"1959-01-01","1979-07-01"),
       x=sub_x(df_us_q,"1959-01-01","1979-07-01")),
  list(label="USA Gran Moderación (trim.)",periodo="Ene 1984 – Oct 2007", freq="Q",
       y=sub_y(df_us_q,"1984-01-01","2007-10-01"),
       x=sub_x(df_us_q,"1984-01-01","2007-10-01")),
  list(label="USA Volcker-Greenspan (trim.)", periodo="1979Q4 – 2005Q4", freq="Q",
       y=sub_y(df_us_q,"1979-10-01","2005-10-01"),
       x=sub_x(df_us_q,"1979-10-01","2005-10-01")),
  list(label="USA Bernanke/Yellen/Powell (trim.)", periodo="2006Q1 – 2025Q2", freq="Q",
       y=sub_y(df_us_q,"2006-01-01","2025-04-01"),
       x=sub_x(df_us_q,"2006-01-01","2025-04-01"))
  
)


# ------------------------------------------------------------
# Calcular
# ------------------------------------------------------------
res <- lapply(muestras, function(m) {
  ri <- persistence_row(m$y)
  ru <- persistence_row(m$x)
  n  <- length(na.omit(m$x))
  list(
    label    = m$label,
    periodo  = m$periodo,
    n        = n,
    i_ar1    = ri$ar1,  i_pval = ri$pval,
    u_ar1    = ru$ar1,  u_pval = ru$pval
  )
})


# ============================================================
# EXPORTAR DIAGNÓSTICOS DE PERSISTENCIA
# ============================================================

persistence_table <- do.call(rbind, lapply(res, function(r) {
  data.frame(
    Sample             = r$label,
    Period             = r$periodo,
    n                  = r$n,
    Inflation_AR1      = r$i_ar1,
    Inflation_ADF_p    = r$i_pval,
    Unemployment_AR1   = r$u_ar1,
    Unemployment_ADF_p = r$u_pval,
    stringsAsFactors   = FALSE
  )
}))

dir.create(
  "Results/Tables/Diagnostics",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  persistence_table,
  "Results/Tables/Diagnostics/persistence_diagnostics.csv",
  row.names = FALSE
)



# ------------------------------------------------------------
# Imprimir tabla limpia
# ------------------------------------------------------------
W <- 97
SEP_D  <- paste0("+", strrep("-", W), "+")
SEP_E  <- paste0("+", strrep("=", W), "+")
SEP_LT <- paste0("|", strrep("-", W), "|")

center <- function(s, w) {
  p <- max(0, w - nchar(s))
  paste0(strrep(" ", floor(p/2)), s, strrep(" ", ceiling(p/2)))
}

cat("\n")
cat(SEP_E, "\n", sep="")
cat(sprintf("| %s |\n", center("Tabla — Persistencia de Inflación y Tasa de Desempleo", W-2)))
cat(SEP_E, "\n", sep="")

cat(sprintf("| %-30s | %-22s | %4s | %13s | %13s |\n",
            "", "", "",
            "  Inflación  ",
            "  Desempleo  "))
cat(sprintf("| %-30s | %-22s | %4s | %6s %6s | %6s %6s |\n",
            "Muestra", "Período", "n",
            "AR(1)", "ADF p",
            "AR(1)", "ADF p"))
cat(SEP_D, "\n", sep="")

for (i in seq_along(res)) {
  r <- res[[i]]
  
  if (i == 5)  cat(SEP_LT, "\n", sep="")  # EA (1-4) → USA mensual (5-9)
  if (i == 10) cat(SEP_D,  "\n", sep="")  # USA mensual → USA trimestral (10-14)
  
  cat(sprintf("| %-30s | %-22s | %4d | %6.3f %6s | %6.3f %6s |\n",
              r$label, r$periodo, r$n,
              r$i_ar1, r$i_pval,
              r$u_ar1, r$u_pval))
}

cat(SEP_D, "\n", sep="")
cat(sprintf("| %-95s |\n",
            "  Notas: ADF con drift. Lag = floor(n^{1/3}). p-valor por tramos MacKinnon."))
cat(sprintf("| %-95s |\n",
            "  AR(1) estimado por MCO: x_t = a + rho * x_{t-1} + e_t"))
cat(SEP_E, "\n\n", sep="")