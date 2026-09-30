# ============================================================
# 01_build_ea_data.R
# INPUT — datos Área Euro (EA-MD-QD)
#
# OUTPUT para IVX (00_ivx_functions.R):
#   df_ea        → data.frame {date, y, x}  longitud T_eff
#   x_full_ea    → vector x completo        longitud T_eff + 1
#
# OUTPUT para XDlasso (06_run_xdlasso.R):
#   data_trans_ea    → panel EA-MD transformado  {date, UNETOT_level, infl_hicp, W...}
#   data_notrans_ea  → panel EA-MD sin transformar (excl. I(2))
#   vars_ea          → nombres controles W (transformado)
#   vars_notrans_ea  → nombres controles W (sin transformar)
#
# NOTA DE ALINEACIÓN:
#   df_ea / x_full_ea  → y ya desfasada (π_{t+1}), para IVX
#   data_trans_ea      → y en nivel (π_t), el lag lo hace main_xdlasso()
# ============================================================

library(readxl)
library(dplyr)

ea_path <- "Data/EA/EA-MD-QD-05-2025/EAdata.xlsx"

# ============================================================
# BLOQUE A: CARGA RAW
# ============================================================

data_raw_ea   <- read_excel(ea_path, sheet = "data") %>% rename(date = Time)
info_ea       <- read_excel(ea_path, sheet = "info")

# Tcodes de EA-MD (columna TR1 en sheet info)
monthly_names_ea <- info_ea %>% filter(Frequency == "M") %>% pull(Name)
tcodes_ea <- info_ea %>%
  filter(Frequency == "M") %>%
  select(Name, TR1) %>%
  tibble::deframe()

cat("=== EA DATA ===\n")
cat("Series mensuales disponibles:", length(monthly_names_ea), "\n")
cat("Verificación UNETOT_EA TR1 =", tcodes_ea["UNETOT_EA"], "(esperado: 4.5)\n")
cat("Verificación HICPOV_EA TR1 =", tcodes_ea["HICPOV_EA"], "(esperado: 2.5)\n")

# ============================================================
# BLOQUE B: OBJETOS PARA IVX
# Limpiar NAs en nivel antes de diff → contrato length(x_full_ea) = nrow(df_ea) + 1
# ============================================================

data_ivx_clean <- data_raw_ea[!is.na(data_raw_ea$HICPOV_EA) &
                               !is.na(data_raw_ea$UNETOT_EA), ]
n_ivx     <- nrow(data_ivx_clean)
x_full_ea <- data_ivx_clean$UNETOT_EA                          # longitud n_ivx = T_eff+1
y_ivx     <- 100 * diff(log(data_ivx_clean$HICPOV_EA))         # longitud n_ivx-1 = T_eff

df_ea <- data.frame(
  date = as.Date(data_ivx_clean$date[1:(n_ivx - 1)]),
  y    = y_ivx,
  x    = x_full_ea[1:(n_ivx - 1)]
)

stopifnot(length(x_full_ea) == nrow(df_ea) + 1)

cat("\n--- IVX objects ---\n")
cat("Rango df_ea:", as.character(min(df_ea$date)),
    "a", as.character(max(df_ea$date)), "\n")
cat("T_eff:", nrow(df_ea), "| x_full_ea:", length(x_full_ea), "\n")

# ============================================================
# BLOQUE C: FUNCIÓN DE TRANSFORMACIÓN EA-MD
# Tcodes EA-MD:
#   0   → nivel
#   2   → 100 * diff(log(x))    (log-diff, produce 1 NA inicial)
#   2.5 → 100 * diff(log(x))    (igual que 2)
#   4   → diff(x)               (primera diferencia, 1 NA inicial)
#   4.5 → diff(x)               (igual que 4)
#   5   → diff(diff(x))         (segunda diferencia, 2 NAs iniciales)
# ============================================================

ea_transform <- function(x, tcode) {
  if (tcode == 0)                  return(x)
  if (tcode == 2 || tcode == 2.5) return(c(NA, 100 * diff(log(x))))
  if (tcode == 4 || tcode == 4.5) return(c(NA, diff(x)))
  if (tcode == 5)                  return(c(NA, NA, diff(diff(x))))
  warning(paste("TCODE desconocido:", tcode))
  return(x)
}

ea_notransform <- function(x, tcode) {
  if (tcode == 0)                  return(x)
  if (tcode == 2 || tcode == 2.5) return(100 * log(x))   # log-nivel
  if (tcode == 4 || tcode == 4.5) return(x)               # nivel directo
  if (tcode == 5)                  return(rep(NA, length(x)))  # I(2) → excluir
  return(x)
}

# ============================================================
# BLOQUE D: PANEL TRANSFORMADO PARA XDLASSO
# Variables excluidas de controles W:
#   UNETOT_EA  → predictor de interés (va separado como d)
#   UNEO25_EA, UNEU25_EA → correlacionadas con predictor
#   HICPOV_EA  → variable dependiente
# ============================================================

vars_excluir_ea <- c("UNETOT_EA", "UNEO25_EA", "UNEU25_EA", "HICPOV_EA")
vars_controles_ea <- setdiff(monthly_names_ea, vars_excluir_ea)

# Transformar controles
controles_trans <- data_raw_ea %>% select(date, all_of(vars_controles_ea))
for (v in vars_controles_ea) {
  controles_trans[[v]] <- ea_transform(data_raw_ea[[v]], tcode = tcodes_ea[v])
}

# y y x para XDlasso (sin desfasar — main_xdlasso hace el lag)
infl_hicp_ea    <- ea_transform(data_raw_ea$HICPOV_EA, tcode = 2.5)
UNETOT_level_ea <- data_raw_ea$UNETOT_EA

# Unir y limpiar
data_trans_ea <- controles_trans %>%
  mutate(
    infl_hicp    = infl_hicp_ea,
    UNETOT_level = UNETOT_level_ea
  ) %>%
  select(date, UNETOT_level, infl_hicp, all_of(vars_controles_ea)) %>%
  select(-any_of("TRNCOG_EA")) %>%   # termina nov-2023, genera NAs al final
  filter(!is.na(infl_hicp)) %>%
  slice(-(1:1)) %>%                   # tcode=5 produce 2 NAs iniciales, slice quita 1 extra
  filter(if_all(everything(), ~ !is.na(.)))

vars_ea <- setdiff(names(data_trans_ea), c("date", "UNETOT_level", "infl_hicp"))

cat("\n--- XDlasso panel (transformado) ---\n")
cat("Dimensiones:", dim(data_trans_ea), "\n")
cat("Fechas:", as.character(data_trans_ea$date[1]),
    "a", as.character(data_trans_ea$date[nrow(data_trans_ea)]), "\n")
cat("Controles W:", length(vars_ea), "\n")
cat("NAs totales:", sum(is.na(data_trans_ea)), "\n")

# ============================================================
# BLOQUE E: PANEL SIN TRANSFORMAR PARA XDLASSO (robustness)
# Variables I(2) (tcode=5) se excluyen en lugar de diferenciarlas
# y sigue siendo 100*diff(log) — necesita ser estacionaria
# ============================================================

vars_i2_ea <- names(tcodes_ea[tcodes_ea == 5])
vars_notrans_excluir <- c(vars_excluir_ea, vars_i2_ea)
vars_notrans_controles <- setdiff(monthly_names_ea, vars_notrans_excluir)

controles_notrans <- data_raw_ea %>% select(date, all_of(vars_notrans_controles))
for (v in vars_notrans_controles) {
  controles_notrans[[v]] <- ea_notransform(data_raw_ea[[v]], tcode = tcodes_ea[v])
}

data_notrans_ea <- controles_notrans %>%
  mutate(
    infl_hicp    = infl_hicp_ea,
    UNETOT_level = UNETOT_level_ea
  ) %>%
  select(date, UNETOT_level, infl_hicp, all_of(vars_notrans_controles)) %>%
  filter(!is.na(infl_hicp)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

vars_notrans_ea <- setdiff(names(data_notrans_ea),
                            c("date", "UNETOT_level", "infl_hicp"))

cat("\n--- XDlasso panel (sin transformar, excl. I(2)) ---\n")
cat("Dimensiones:", dim(data_notrans_ea), "\n")
cat("Fechas:", as.character(data_notrans_ea$date[1]),
    "a", as.character(data_notrans_ea$date[nrow(data_notrans_ea)]), "\n")
cat("Controles W:", length(vars_notrans_ea), "\n")
cat("Variables I(2) excluidas:", paste(vars_i2_ea, collapse = ", "), "\n")
cat("NAs totales:", sum(is.na(data_notrans_ea)), "\n\n")

cat("✅ 01_build_ea_data.R OK\n")
cat("   Objetos IVX:     df_ea, x_full_ea\n")
cat("   Objetos XDlasso: data_trans_ea, vars_ea\n")
cat("                    data_notrans_ea, vars_notrans_ea\n")
