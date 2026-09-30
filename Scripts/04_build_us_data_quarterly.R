# ============================================================
# 02b_build_us_data_quarterly.R
#
# PROPÓSITO
# ---------
# Construir un panel TRIMESTRAL para EE.UU. a partir de FRED-MD (mensual),
# en paralelismo exacto con 01b_build_ea_data_quarterly.R para el Área Euro.
#
# Por qué hacemos esto:
#   02_build_us_data.R produce datos mensuales (T ≈ 800).
#   Al agregar a trimestral obtenemos T ≈ 260 observaciones trimestrales,
#   con el mismo universo de ~110 controles W, permitiendo:
#     (1) Comparación directa mensual vs. trimestral en EE.UU.
#     (2) Comparación simétrica con el panel trimestral del Área Euro
#     (3) Robustez: los resultados IVX/XDlasso no deberían cambiar
#         drásticamente al cambiar frecuencia
#
# INPUT
# -----
#   FRED-MD vintage CSV (us_path)
#   Fila 1: nombres de columnas
#   Fila 2: tcodes ("Transform:", 5, 5, ...)
#
# OUTPUT para IVX (00_ivx_functions.R):
#   df_us_q        → data.frame {date, y, x}  longitud T_eff
#   x_full_us_q    → vector x completo        longitud T_eff + 1
#
# OUTPUT para XDlasso (06_run_xdlasso.R):
#   data_trans_us_q    → panel trimestral transformado {date, UNRATE_level, infl_cpi, W...}
#   data_notrans_us_q  → panel trimestral sin transformar (excl. I(2))
#   vars_us_q          → nombres controles W (transformado)
#   vars_notrans_us_q  → nombres controles W (sin transformar)
#
# NOTA DE ALINEACIÓN (igual que versión mensual):
#   df_us_q / x_full_us_q  → y ya desfasada (π_{t+1}), para IVX
#   data_trans_us_q        → y en nivel (π_t), el lag lo hace main_xdlasso()
#
# ============================================================
# DECISIONES DE IMPLEMENTACIÓN
# ============================================================
#
# [D1] AGREGACIÓN DE SERIES MENSUALES → TRIMESTRAL
#
#   FRED-MD es mensual. Para pasar a trimestral seguimos el mismo criterio
#   que en 01b para EA-MD, usando el principio general:
#     → Todos los controles: MEDIA de los 3 meses del trimestre sobre el NIVEL
#       crudo (antes de transformar). Si algún mes es NA, el trimestre es NA.
#     → UNRATE (tcode=2, nivel): media de los 3 meses → luego diff() trimestral
#     → CPIAUCSL (tcode=6, índice de precios): media del índice crudo del
#       trimestre → luego 100*diff(diff(log())) trimestral
#
#   IMPORTANTE: NO usamos fbi::fredmd(transform=TRUE) para el panel trimestral
#   porque fbi aplica las transformaciones en frecuencia mensual. Aquí
#   necesitamos agregar PRIMERO y transformar DESPUÉS (igual que EA quarterly).
#   Por eso replicamos los tcodes FRED manualmente (función us_transform_q).
#
# [D2] INFLACIÓN TRIMESTRAL
#
#   CPIAUCSL tiene tcode=6 → diff(diff(log(x))) en mensual.
#   En trimestral usamos: π_Q(t) = 100 * diff(log(CPI_Q(t)))
#   donde CPI_Q(t) = media del índice crudo en el trimestre t.
#   Esto equivale a la tasa de crecimiento log del nivel de precios trimestral,
#   análogo a la convención EA quarterly con HICPOV_EA.
#
#   Nota: el tcode original 6 implica segunda diferencia del log (aceleración
#   de la inflación). Para la variable dependiente y usamos directamente la
#   primera diferencia del log del índice trimestral, que es la inflación
#   trimestral estándar. Esto es coherente con el uso en el script mensual
#   (bloque B usa diff(log(CPI)) a pesar de tcode=6).
#
# [D3] TCODES FRED-MD vs EA-MD
#
#   FRED-MD usa tcodes 1-7:
#     1 → nivel                         (sin transformar)
#     2 → diff(x)                       (primera diferencia)
#     3 → diff(diff(x))                 (segunda diferencia)  ← I(2)
#     4 → log(x)                        (log-nivel)
#     5 → diff(log(x))                  (log-diferencia)
#     6 → diff(diff(log(x)))            (segunda diferencia del log) ← I(2)
#     7 → diff(x/lag(x) - 1)            (cambio porcentual de la tasa de cambio)
#
#   Los tcodes se aplican sobre series ya agregadas a trimestral.
#   Las I(2) (tcodes 3 y 6) se excluyen en el panel sin transformar (bloque E).
#
# [D4] ALINEACIÓN DE FECHAS
#
#   Al agregar mensual a trimestral, cada grupo de 3 meses recibe
#   la fecha del PRIMER mes del trimestre via floor_date(..., "quarter"):
#     Ene/Feb/Mar → YYYY-01-01
#     Abr/May/Jun → YYYY-04-01
#     Jul/Ago/Sep → YYYY-07-01
#     Oct/Nov/Dic → YYYY-10-01
#   Esto es consistente con la convención de 01b_build_ea_data_quarterly.R.
#
# ============================================================

library(dplyr)
library(lubridate)

us_path <- "Data/US/2026-04-MD.csv"

# ============================================================
# BLOQUE A: CARGA RAW
#
# Misma lógica que 02_build_us_data.R:
#   Fila "Transform:" contiene los tcodes
#   Resto de filas son los datos mensuales
# ============================================================

raw_all      <- read.csv(us_path, header = TRUE, stringsAsFactors = FALSE)
tcode_row    <- raw_all[raw_all$sasdate == "Transform:", ]
tcodes_us    <- as.numeric(tcode_row[1, -1])
var_names_us_q      <- names(raw_all)[-1]
var_names_us_q_norm <- make.names(var_names_us_q)  # normaliza S&P 500 → S.P.500, etc.
names(tcodes_us)    <- var_names_us_q_norm         # tcodes indexados por nombre normalizado

data_raw_us <- raw_all[raw_all$sasdate != "Transform:", ]
data_raw_us$date <- as.Date(data_raw_us$sasdate, format = "%m/%d/%Y")
for (v in var_names_us_q) {
  data_raw_us[[v]] <- as.numeric(data_raw_us[[v]])
}
# Renombrar columnas con nombres normalizados (una sola vez, aquí)
names(data_raw_us)[match(var_names_us_q, names(data_raw_us))] <- var_names_us_q_norm

cat("=== US DATA TRIMESTRAL (FRED-MD) ===\n")
cat("Series disponibles:", length(var_names_us_q), "\n")
cat("Verificación CPIAUCSL tcode =", tcodes_us["CPIAUCSL"], "(esperado: 6)\n")
cat("Verificación UNRATE   tcode =", tcodes_us["UNRATE"],   "(esperado: 2)\n")

# ============================================================
# BLOQUE B: FUNCIONES DE TRANSFORMACIÓN
#
# us_transform_q() aplica los tcodes FRED sobre series ya agregadas
# a frecuencia trimestral. Paralelo exacto a ea_transform_q().
#
# us_notransform_q() devuelve nivel o log-nivel para el panel de
# robustez (bloque E). I(2) devuelven NA → excluidas via filter().
#
# NOTA sobre tcode=7 (diff(x/lag(x)-1)):
#   En frecuencia mensual, tcode=7 = diff de la tasa de crecimiento mensual.
#   Al agregar a trimestral, primero promediamos el nivel crudo y luego
#   aplicamos la misma fórmula sobre el nivel trimestral.
#   Es una aproximación razonable y coherente con el tratamiento del resto.
# ============================================================

us_transform_q <- function(x, tcode) {
  if (tcode == 1) return(x)                                   # nivel
  if (tcode == 2) return(c(NA, diff(x)))                      # primera diferencia
  if (tcode == 3) return(c(NA, NA, diff(diff(x))))            # segunda diferencia (I(2))
  if (tcode == 4) return(log(x))                              # log-nivel (sin diff)
  if (tcode == 5) return(c(NA, diff(log(x))))                 # log-diferencia
  if (tcode == 6) return(c(NA, NA, diff(diff(log(x)))))       # segunda diff del log (I(2))
  if (tcode == 7) {                                           # diff(x/lag(x) - 1)
    growth <- x / c(NA, x[-length(x)]) - 1
    return(c(NA, diff(growth)))
  }
  warning(paste("TCODE desconocido:", tcode))
  return(x)
}

us_notransform_q <- function(x, tcode) {
  if (tcode == 1) return(x)                          # nivel directo
  if (tcode == 2) return(x)                          # nivel directo
  if (tcode == 3) return(rep(NA, length(x)))         # I(2) → excluir
  if (tcode == 4) return(log(x))                     # log-nivel
  if (tcode == 5) return(log(x))                     # log-nivel
  if (tcode == 6) return(rep(NA, length(x)))         # I(2) → excluir
  if (tcode == 7) return(x)                          # nivel
  return(x)
}

# ============================================================
# BLOQUE C: AGREGACIÓN DE SERIES MENSUALES → TRIMESTRAL
#
# Para CADA serie mensual de FRED-MD:
#   1. Tomamos el nivel crudo (sin transformar)
#   2. Agrupamos por trimestre (floor_date a "quarter")
#   3. Calculamos la MEDIA de los 3 meses
#      (na.rm = FALSE: si algún mes es NA, el trimestre es NA)
#
# Resultado: data_agg_q con una fila por trimestre y todas las
# series en nivel crudo, listas para transformar en bloques D y E.
# ============================================================

cat("\n--- Agregando series mensuales a trimestral ---\n")

data_monthly <- data_raw_us %>%
  mutate(quarter = floor_date(date, "quarter"))

# Calcular media trimestral para todas las variables numéricas
data_agg_q <- data_monthly %>%
  group_by(quarter) %>%
  summarise(
    across(all_of(var_names_us_q_norm), ~ if (sum(!is.na(.)) == 3) mean(., na.rm = TRUE) else NA_real_),
    .groups = "drop"
  ) %>%
  rename(date = quarter)

# Renombrar columnas de data_agg_q con nombres normalizados
# (make.names aplicado sobre var_names_us_q_norm que ya está definido en bloque A)
names(data_agg_q)[match(var_names_us_q, names(data_agg_q))] <- var_names_us_q_norm

cat("Dimensiones panel agregado:", dim(data_agg_q), "\n")
cat("Rango:", as.character(min(data_agg_q$date)),
    "a", as.character(max(data_agg_q$date)), "\n")

# ============================================================
# BLOQUE D: OBJETOS PARA IVX — VERSIÓN TRIMESTRAL
#
# y = inflación trimestral CPI = 100 * diff(log(CPIAUCSL_Q))
#     donde CPIAUCSL_Q es la media del índice crudo del trimestre
# x = UNRATE en nivel = media de los 3 meses del trimestre
#
# Alineación temporal (igual que script mensual y EA quarterly):
#   x[t] predice y[t+1]
#   x_full_us_q tiene longitud T_eff + 1
#   df_us_q contiene los pares (y[t+1], x[t]) para t = 1, ..., T_eff
# ============================================================

cat("\n--- Construyendo objetos IVX trimestrales ---\n")

data_ivx_q <- data_agg_q %>%
  filter(!is.na(CPIAUCSL) & !is.na(UNRATE))

n_ivx_q     <- nrow(data_ivx_q)
x_full_us_q <- data_ivx_q$UNRATE                          # longitud = T_eff + 1
y_ivx_q     <- 100 * diff(log(data_ivx_q$CPIAUCSL))       # longitud = T_eff

df_us_q <- data.frame(
  date = data_ivx_q$date[1:(n_ivx_q - 1)],
  y    = y_ivx_q,
  x    = x_full_us_q[1:(n_ivx_q - 1)]
)

stopifnot(length(x_full_us_q) == nrow(df_us_q) + 1)

cat("Rango df_us_q:", as.character(min(df_us_q$date)),
    "a", as.character(max(df_us_q$date)), "\n")
cat("T_eff:", nrow(df_us_q), "| x_full_us_q:", length(x_full_us_q), "\n")

# ============================================================
# BLOQUE E: PANEL TRANSFORMADO PARA XDLASSO — VERSIÓN TRIMESTRAL
#
# Construimos la matriz de controles W con todas las series de
# FRED-MD agregadas a trimestral, transformadas a estacionarias
# según sus tcodes FRED (aplicados sobre el nivel trimestral).
#
# Variables excluidas de W (misma lógica que 02_build_us_data.R):
#   UNRATE    → predictor de interés (d), va separado
#   CPIAUCSL  → variable dependiente (y)
#
# Selección final de controles: solo las series sin NAs en el
# rango de muestra, usando select_if(~ !any(is.na(.))) igual
# que en el script mensual de Gao et al.
#
# RANGO: filtramos ANTES de select_if para no penalizar series
# que empiezan tarde (NAs iniciales) o trimestres incompletos al
# final (Q4-2025 y Q1-2026 solo tienen 2 meses en FRED-MD 2026-03).
#   date_start_q = "1960-01-01"  → primer trimestre completo con
#                                   historia suficiente en FRED-MD
#   date_end_q   = "2025-07-01"  → último trimestre con 3 meses OK
#                                   (verificado: Q4-2025 solo tiene 2)
#
# El panel NO tiene desfase — main_xdlasso() hace el lag internamente.
# ============================================================

date_start_q <- as.Date("1960-01-01")
date_end_q   <- as.Date("2025-07-01")

cat("\n--- Construyendo panel XDlasso trimestral (transformado) ---\n")

vars_excluir_q_us <- c("UNRATE", "CPIAUCSL")
# Restringir al mismo universo de controles que sobrevive en el panel mensual
# (vars_us definido en 02_build_us_data.R). Garantiza comparabilidad directa
# mensual vs. trimestral: mismas 110 series, distinta frecuencia.
if (!exists("vars_us")) stop("vars_us no encontrado — ejecuta primero 02_build_us_data.R")
vars_controles_q_us <- intersect(var_names_us_q_norm, make.names(vars_us))

# Calcular y e x para XDlasso (sin desfasar — main_xdlasso hace el lag)
infl_cpi_q     <- c(NA, 100 * diff(log(data_agg_q$CPIAUCSL)))
UNRATE_level_q <- data_agg_q$UNRATE

# Transformar cada control según su tcode FRED, aplicado sobre nivel trimestral
controles_trans_q <- data_agg_q %>%
  select(date, all_of(vars_controles_q_us))

for (v in vars_controles_q_us) {
  controles_trans_q[[v]] <- us_transform_q(data_agg_q[[v]], tcode = tcodes_us[v])
}

# Unir y, x y controles; limpiar NAs
data_trans_us_q <- controles_trans_q %>%
  mutate(
    infl_cpi     = infl_cpi_q,
    UNRATE_level = UNRATE_level_q
  ) %>%
  select(date, UNRATE_level, infl_cpi, all_of(vars_controles_q_us)) %>%
  filter(!is.na(infl_cpi)) %>%                         # 1 NA inicial por diff(log)
  filter(date >= date_start_q, date <= date_end_q) %>% # rango con historia completa
  select_if(~ !any(is.na(.)))                          # excluir series con cualquier NA

vars_us_q <- setdiff(
  names(data_trans_us_q),
  c("date", "UNRATE_level", "infl_cpi")
)

cat("Dimensiones:", dim(data_trans_us_q), "\n")
cat("Fechas:", as.character(data_trans_us_q$date[1]),
    "a", as.character(data_trans_us_q$date[nrow(data_trans_us_q)]), "\n")
cat("Controles W:", length(vars_us_q), "\n")
cat("NAs totales:", sum(is.na(data_trans_us_q)), "\n")

# ============================================================
# BLOQUE F: PANEL SIN TRANSFORMAR PARA XDLASSO (robustez)
#
# Igual que bloque E del script mensual y bloque H de EA quarterly:
# series en nivel o log-nivel, excluyendo I(2) (tcodes 3 y 6).
#
# Este panel sirve como ejercicio de robustez: si XDlasso cambia
# al usar niveles vs. diferencias, la inferencia es sensible a
# la especificación de los controles.
# ============================================================

cat("\n--- Construyendo panel XDlasso trimestral (sin transformar) ---\n")

vars_i2_us <- names(tcodes_us[tcodes_us %in% c(3, 6)])

# Misma restricción: solo controles que sobreviven en el panel mensual sin transformar
# (vars_notrans_us definido en 02_build_us_data.R)
if (!exists("vars_notrans_us")) stop("vars_notrans_us no encontrado — ejecuta primero 02_build_us_data.R")
vars_notrans_excluir_q <- c(vars_excluir_q_us, vars_i2_us)
vars_notrans_controles_q <- intersect(var_names_us_q_norm, make.names(vars_notrans_us))

controles_notrans_q <- data_agg_q %>%
  select(date, all_of(vars_notrans_controles_q))

for (v in vars_notrans_controles_q) {
  controles_notrans_q[[v]] <- us_notransform_q(data_agg_q[[v]], tcode = tcodes_us[v])
}

data_notrans_us_q <- controles_notrans_q %>%
  mutate(
    infl_cpi     = infl_cpi_q,
    UNRATE_level = UNRATE_level_q
  ) %>%
  select(date, UNRATE_level, infl_cpi, all_of(vars_notrans_controles_q)) %>%
  filter(!is.na(infl_cpi)) %>%
  filter(date >= date_start_q, date <= date_end_q) %>% # mismo rango que bloque E
  select_if(~ !any(is.na(.)))

vars_notrans_us_q <- setdiff(
  names(data_notrans_us_q),
  c("date", "UNRATE_level", "infl_cpi")
)

cat("Dimensiones:", dim(data_notrans_us_q), "\n")
cat("Fechas:", as.character(data_notrans_us_q$date[1]),
    "a", as.character(data_notrans_us_q$date[nrow(data_notrans_us_q)]), "\n")
cat("Controles W:", length(vars_notrans_us_q), "\n")
cat("Variables I(2) excluidas:", paste(vars_i2_us, collapse = ", "), "\n")
cat("NAs totales:", sum(is.na(data_notrans_us_q)), "\n\n")

cat("✅ 02b_build_us_data_quarterly.R OK\n")
cat("   Objetos IVX:     df_us_q, x_full_us_q\n")
cat("   Objetos XDlasso: data_trans_us_q, vars_us_q\n")
cat("                    data_notrans_us_q, vars_notrans_us_q\n")
