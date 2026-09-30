# ============================================================
# 02_build_us_data.R
# INPUT  — FRED-MD vintage (us_path)
#
# OUTPUT para IVX (00_ivx_functions.R):
#   df_us        → data.frame {date, y, x}  longitud T_eff
#   x_full_us    → vector x completo        longitud T_eff + 1
#
# OUTPUT para XDlasso (06_run_xdlasso.R):
#   data_trans_us    → panel transformado   {date, UNRATE_level, infl_cpi, W...}
#   data_notrans_us  → panel sin transformar, excl. I(2)
#   vars_us          → nombres controles W (transformado)
#   vars_notrans_us  → nombres controles W (sin transformar)
#
# NOTA DE ALINEACIÓN:
#   df_us / x_full_us  → y ya desfasada (π_{t+1}), para IVX
#   data_trans_us      → y en nivel (π_t), el lag lo hace main_xdlasso()
#
# LÓGICA DE CONTROLES W (bloque D):
#   Replica exactamente el repo Gao et al. (github.com/zhan-gao/lasso_inference):
#   fbi::fredmd(transform=TRUE) → intersect(varlist$fred) → select_if(!any(is.na(.)))
#   Con vintage marzo 2026, end_date sept 2025 → 110 controles W
#   (Gao usa vintage agosto 2025, end_date abril 2025 → 110 controles también)
#
# LÓGICA BLOQUE E (sin transformar):
#   Excluye variables I(2) por tcode (3 y 6) y variables sin historia pre-1960.
#   Usa la misma lista de exclusión que el bloque D para UNRATE y CPIAUCSL.
#
# TCODES FRED-MD:
#   1 → nivel
#   2 → diff(x)
#   3 → diff(diff(x))         ← I(2), excluido en panel sin transformar
#   4 → log(x)
#   5 → diff(log(x))
#   6 → diff(diff(log(x)))    ← I(2), excluido en panel sin transformar
#   7 → diff(x/lag(x) - 1)
# ============================================================

library(dplyr)
library(fbi)

us_path <- "Data/US/2026-04-MD.csv"

# ============================================================
# BLOQUE A: CARGA RAW
# Fila 1: nombres | Fila 2: tcodes (Transform:, 5, 5, ...)
# Necesario para IVX (bloque B) y panel sin transformar (bloque E)
# ============================================================

raw_all      <- read.csv(us_path, header = TRUE, stringsAsFactors = FALSE)
tcode_row    <- raw_all[raw_all$sasdate == "Transform:", ]
tcodes_us    <- as.numeric(tcode_row[1, -1])
var_names_us <- names(raw_all)[-1]
names(tcodes_us) <- var_names_us

data_raw_us <- raw_all[raw_all$sasdate != "Transform:", ]
data_raw_us$date <- as.Date(data_raw_us$sasdate, format = "%m/%d/%Y")
for (v in var_names_us) {
  data_raw_us[[v]] <- as.numeric(data_raw_us[[v]])
}

cat("=== US DATA (FRED-MD) ===\n")
cat("Series disponibles:", length(var_names_us), "\n")
cat("Verificación CPIAUCSL tcode =", tcodes_us["CPIAUCSL"], "(esperado: 6)\n")
cat("Verificación UNRATE   tcode =", tcodes_us["UNRATE"],   "(esperado: 2)\n")

# ============================================================
# BLOQUE B: OBJETOS PARA IVX
# y = 100 * diff(log(CPI))  → inflación mensual en %
# x = UNRATE en nivel       → desempleo sin transformar
# Alineación: x[t] predice y[t+1], por eso x_full tiene longitud T_eff + 1
# ============================================================

data_ivx_clean <- data_raw_us[!is.na(data_raw_us$CPIAUCSL) &
                               !is.na(data_raw_us$UNRATE), ]
n_ivx     <- nrow(data_ivx_clean)
x_full_us <- data_ivx_clean$UNRATE
y_ivx     <- 100 * diff(log(data_ivx_clean$CPIAUCSL))

df_us <- data.frame(
  date = data_ivx_clean$date[1:(n_ivx - 1)],
  y    = y_ivx,
  x    = x_full_us[1:(n_ivx - 1)]
)

stopifnot(length(x_full_us) == nrow(df_us) + 1)

cat("\n--- IVX objects ---\n")
cat("Rango df_us:", as.character(min(df_us$date)),
    "a", as.character(max(df_us$date)), "\n")
cat("T_eff:", nrow(df_us), "| x_full_us:", length(x_full_us), "\n")

# ============================================================
# BLOQUE C: FUNCIÓN DE TRANSFORMACIÓN MANUAL (solo para bloque E)
# No se usa en bloque D — fbi::fredmd aplica tcodes internamente
# ============================================================

us_notransform <- function(x, tcode) {
  # Devuelve la serie en nivel apropiado para el panel sin transformar.
  # I(2) (tcodes 3 y 6) → NA para exclusión posterior via filter().
  if (tcode == 1) return(x)
  if (tcode == 2) return(x)                     # nivel directo
  if (tcode == 3) return(rep(NA, length(x)))    # I(2) → excluir
  if (tcode == 4) return(log(x))                # log-nivel
  if (tcode == 5) return(log(x))                # log-nivel
  if (tcode == 6) return(rep(NA, length(x)))    # I(2) → excluir
  if (tcode == 7) return(x)                     # nivel
  return(x)
}

# ============================================================
# BLOQUE D
# ============================================================

varlist_fred <- fbi::fredmd_description

data_fbi_trans <- fbi::fredmd(
  us_path,
  date_start = as.Date("1959-11-01"),
  date_end   = as.Date("2025-09-01"),
  transform  = TRUE
)

data_fbi_raw <- fbi::fredmd(
  us_path,
  date_start = as.Date("1959-11-01"),
  date_end   = as.Date("2025-09-01"),
  transform  = FALSE
)

infl_cpi_fbi     <- c(NA, diff(log(data_fbi_raw$CPIAUCSL)) * 100)
UNRATE_level_fbi <- data_fbi_raw$UNRATE

vars_fbi <- intersect(colnames(data_fbi_trans), varlist_fred$fred)

data_trans_us <- data_fbi_trans %>%
  as_tibble() %>%
  select(all_of(c("date", vars_fbi))) %>%
  select_if(~ !any(is.na(.))) %>%
  mutate(
    infl_cpi     = infl_cpi_fbi[match(.$date, data_fbi_raw$date)],
    UNRATE_level = UNRATE_level_fbi[match(.$date, data_fbi_raw$date)]
  ) %>%
  filter(!is.na(infl_cpi)) %>%
  select(date, UNRATE_level, infl_cpi, everything())

vars_us <- setdiff(
  names(data_trans_us),
  c("date", "UNRATE_level", "infl_cpi", "UNRATE", "CPIAUCSL")
)

controls_110      <- vars_us
controls_110_norm <- make.names(controls_110)  # normaliza S&P 500 → S.P.500, etc.

cat("\n--- XDlasso panel (transformado) ---\n")
cat("Dimensiones:", dim(data_trans_us), "\n")
cat("Fechas:", as.character(data_trans_us$date[1]),
    "a", as.character(data_trans_us$date[nrow(data_trans_us)]), "\n")
cat("Controles W:", length(vars_us), "\n")
cat("NAs totales:", sum(is.na(data_trans_us)), "\n")

# ============================================================
# BLOQUE E
# ============================================================
# ============================================================
# BLOQUE E — Sin transformar, excl. I(2) y excl. NAs
# ============================================================

vars_i2_us      <- names(tcodes_us[tcodes_us %in% c(3, 6)])
vars_i2_us_norm <- make.names(vars_i2_us)

# Universo completo - I(2) - UNRATE - CPIAUCSL
vars_candidatos_notrans <- setdiff(
  make.names(var_names_us),
  make.names(c(vars_i2_us, "CPIAUCSL", "UNRATE"))
)

# Construir panel con todas las candidatas
controles_notrans_us <- data.frame(date = as.Date(data_raw_us$sasdate, 
                                                  format = "%m/%d/%Y"))
controles_notrans_us <- data_raw_us %>%
  select(date) %>%
  mutate(date = as.Date(date))

for (v in vars_candidatos_notrans) {
  v_orig <- var_names_us[make.names(var_names_us) == v]
  controles_notrans_us[[v]] <- us_notransform(
    data_raw_us[[v_orig]], tcode = tcodes_us[v_orig]
  )
}

infl_cpi_notrans     <- c(NA, diff(log(data_raw_us$CPIAUCSL)) * 100)
UNRATE_level_notrans <- data_raw_us$UNRATE

data_notrans_us <- controles_notrans_us %>%
  mutate(
    infl_cpi     = infl_cpi_notrans,
    UNRATE_level = UNRATE_level_notrans,
    date         = as.Date(date)
  ) %>%
  select(date, UNRATE_level, infl_cpi, all_of(vars_candidatos_notrans)) %>%
  filter(date >= as.Date("1960-01-01"),
         date <= as.Date("2025-04-01")) %>%
  filter(!is.na(infl_cpi)) %>%
  select_if(~ !any(is.na(.)))   # ← igual que bloque D

vars_notrans_us <- setdiff(
  names(data_notrans_us),
  c("date", "UNRATE_level", "infl_cpi")
)

cat("\n--- XDlasso panel (sin transformar, excl. I(2)) ---\n")
cat("Dimensiones:", dim(data_notrans_us), "\n")
cat("Fechas:", as.character(data_notrans_us$date[1]),
    "a", as.character(data_notrans_us$date[nrow(data_notrans_us)]), "\n")
cat("Controles W:", length(vars_notrans_us), "\n")
cat("NAs totales:", sum(is.na(data_notrans_us)), "\n\n")# ============================================================
# BLOQUE E — Sin transformar, excl. I(2) y excl. NAs
# ============================================================

vars_i2_us      <- names(tcodes_us[tcodes_us %in% c(3, 6)])
vars_i2_us_norm <- make.names(vars_i2_us)

# Universo completo - I(2) - UNRATE - CPIAUCSL
vars_candidatos_notrans <- setdiff(
  make.names(var_names_us),
  make.names(c(vars_i2_us, "CPIAUCSL", "UNRATE"))
)

# Construir panel con todas las candidatas
controles_notrans_us <- data.frame(date = as.Date(data_raw_us$sasdate, 
                                                  format = "%m/%d/%Y"))
controles_notrans_us <- data_raw_us %>%
  select(date) %>%
  mutate(date = as.Date(date))

for (v in vars_candidatos_notrans) {
  v_orig <- var_names_us[make.names(var_names_us) == v]
  controles_notrans_us[[v]] <- us_notransform(
    data_raw_us[[v_orig]], tcode = tcodes_us[v_orig]
  )
}

infl_cpi_notrans     <- c(NA, diff(log(data_raw_us$CPIAUCSL)) * 100)
UNRATE_level_notrans <- data_raw_us$UNRATE

data_notrans_us <- controles_notrans_us %>%
  mutate(
    infl_cpi     = infl_cpi_notrans,
    UNRATE_level = UNRATE_level_notrans,
    date         = as.Date(date)
  ) %>%
  select(date, UNRATE_level, infl_cpi, all_of(vars_candidatos_notrans)) %>%
  filter(date >= as.Date("1960-01-01"),
         date <= as.Date("2025-04-01")) %>%
  filter(!is.na(infl_cpi)) %>%
  select_if(~ !any(is.na(.)))   # ← igual que bloque D

vars_notrans_us <- setdiff(
  names(data_notrans_us),
  c("date", "UNRATE_level", "infl_cpi")
)

cat("\n--- XDlasso panel (sin transformar, excl. I(2)) ---\n")
cat("Dimensiones:", dim(data_notrans_us), "\n")
cat("Fechas:", as.character(data_notrans_us$date[1]),
    "a", as.character(data_notrans_us$date[nrow(data_notrans_us)]), "\n")
cat("Controles W:", length(vars_notrans_us), "\n")
cat("NAs totales:", sum(is.na(data_notrans_us)), "\n\n")