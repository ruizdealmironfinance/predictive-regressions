# ============================================================
# 01b_build_ea_data_quarterly.R
#
# PROPÓSITO
# ---------
# Construir un panel TRIMESTRAL combinando:
#   (A) Series mensuales del EA-MD-QD → agregadas a trimestral
#   (B) Series trimestrales nativas del EA-MD-QD → usadas directamente
#
# Por qué hacemos esto:
#   El script mensual (01_build_ea_data.R) produjo solo ~42 controles W,
#   insuficientes para justificar el uso de Double Lasso (high-dimensionality
#   requiere p >> n). Al pasar a frecuencia trimestral:
#     - n ≈ 100 observaciones (1999 Q1 → 2025 Q1)
#     - p ≈ 100+ controles W (mensuales agregadas + trimestrales nativas)
#   Con p > n el Lasso está justificado metodológicamente.
#
# INPUT
# -----
#   EAdata.xlsx  →  sheet "data": todas las series (mensuales y trimestrales)
#                   sheet "info": metadatos (nombre, frecuencia, tcodes)
#
# OUTPUT para IVX (00_ivx_functions.R):
#   df_ea_q      → data.frame {date, y, x}   longitud T_eff
#   x_full_ea_q  → vector x completo         longitud T_eff + 1
#
# OUTPUT para XDlasso (06_run_xdlasso.R):
#   data_trans_ea_q    → panel trimestral transformado {date, UNETOT_level, infl_hicp, W...}
#   data_notrans_ea_q  → panel trimestral sin transformar (excl. I(2))
#   vars_ea_q          → nombres controles W (transformado)
#   vars_notrans_ea_q  → nombres controles W (sin transformar)
#
# NOTA DE ALINEACIÓN (igual que versión mensual):
#   df_ea_q / x_full_ea_q  → y ya desfasada (π_{t+1}), para IVX
#   data_trans_ea_q        → y en nivel (π_t), el lag lo hace main_xdlasso()
#
# ============================================================
# DECISIONES DE IMPLEMENTACIÓN
# ============================================================
#
# [D1] AGREGACIÓN DE SERIES MENSUALES → TRIMESTRAL
#
#   El Excel EA-MD-QD tiene UNA FILA POR MES para todas las series.
#   Las series trimestrales nativas solo tienen valor en el último mes
#   del trimestre (marzo, junio, septiembre, diciembre) y NA en los demás.
#   Las series mensuales tienen valor en todos los meses.
#
#   Para pasar series mensuales a trimestral, seguimos el criterio oficial
#   del ReadME del dataset (opción QM):
#     → FLUJOS (tasas, índices, tipos de interés): MEDIA de los 3 meses
#     → STOCKS (agregados monetarios, niveles financieros): SUMA de los 3 meses
#
#   Para UNETOT (tasa de desempleo, flujo) → media
#   Para HICPOV (índice de precios, nivel) → media del índice,
#     luego inflación trimestral = 100 * diff(log(HICPOV_trimestral))
#     [NO se promedian inflaciones mensuales — se calcula sobre el índice agregado]
#
#   CONDICIÓN: los 3 meses del trimestre deben estar disponibles.
#   Si algún mes es NA → el trimestre entero es NA (na.rm = FALSE).
#   Esto es consistente con el criterio del ReadME:
#   "provided all months related to a quarter are available as data points"
#
# [D2] INFLACIÓN TRIMESTRAL
#
#   Calculamos: π_Q(t) = 100 * (log(HICPOV_Q(t)) - log(HICPOV_Q(t-1)))
#   donde HICPOV_Q(t) es la media del índice en el trimestre t.
#
#   Alternativa descartada: promediar inflaciones mensuales.
#   Razón: la convención estándar en macroeconomía (BCE, Eurostat) mide
#   el cambio del nivel de precios promedio del trimestre respecto al anterior.
#   Promediar tasas de cambio introduce una distorsión sutil y no es estándar.
#
# [D3] ALINEACIÓN DE FECHAS PARA EL JOIN
#
#   Al agregar mensuales a trimestral, cada grupo de 3 meses recibe
#   la fecha del PRIMER mes del trimestre usando floor_date(..., "quarter"):
#     Ene/Feb/Mar 2000 → 2000-01-01
#     Abr/May/Jun 2000 → 2000-04-01
#     ...
#   Las series trimestrales nativas (con datos en mar/jun/sep/dic) también
#   se re-etiquetan con floor_date para que las fechas coincidan en el join.
#   Esto garantiza que ambas fuentes queden alineadas temporalmente.
#
# [D4] SERIES I(2)
#
#   Misma lógica que en el script mensual:
#   - Panel transformado (D): series I(2) con tcode=5 se diferencian dos veces
#   - Panel sin transformar (E): series I(2) se excluyen directamente
#     porque no tienen representación en nivel estacionaria válida.
#
# ============================================================

library(readxl)
library(dplyr)
library(lubridate)

ea_path <- "Data/EA/EA-MD-QD-05-2025/EAdata.xlsx"

# ============================================================
# BLOQUE A: CARGA RAW
#
# Leemos las dos hojas del Excel:
#   "data" → todas las series (filas = meses, columnas = variables)
#   "info" → metadatos de cada variable (nombre, frecuencia, tcode)
#
# IMPORTANTE: el Excel tiene una fila por mes para TODAS las series,
# incluyendo las trimestrales. Las trimestrales tienen NAs en los
# meses que no son fin de trimestre (solo tienen dato en mar/jun/sep/dic).
# ============================================================

data_raw_ea <- read_excel(ea_path, sheet = "data") %>%
  rename(date = Time) %>%
  mutate(date = as.Date(date))

info_ea <- read_excel(ea_path, sheet = "info")

# Separar metadatos por frecuencia
monthly_names_ea   <- info_ea %>% filter(Frequency == "M") %>% pull(Name)
quarterly_names_ea <- info_ea %>% filter(Frequency == "Q") %>% pull(Name)

# Tcodes para series mensuales (usamos TR1 = light transformation, igual que script mensual)
tcodes_monthly <- info_ea %>%
  filter(Frequency == "M") %>%
  select(Name, TR1) %>%
  tibble::deframe()

# Tcodes para series trimestrales
tcodes_quarterly <- info_ea %>%
  filter(Frequency == "Q") %>%
  select(Name, TR1) %>%
  tibble::deframe()

cat("=== EA DATA TRIMESTRAL ===\n")
cat("Series mensuales en info:", length(monthly_names_ea), "\n")
cat("Series trimestrales en info:", length(quarterly_names_ea), "\n")
cat("Verificación UNETOT_EA tcode =", tcodes_monthly["UNETOT_EA"],
    "(esperado: 0 — nivel)\n")
cat("Verificación HICPOV_EA tcode =", tcodes_monthly["HICPOV_EA"],
    "(esperado: 2.5 — log-diff)\n")

# ============================================================
# BLOQUE B: FUNCIONES DE TRANSFORMACIÓN
#
# Necesitamos dos funciones, igual que en el script mensual:
#
# ea_transform_q()    → aplica transformación completa (para panel D)
#   Produce series estacionarias según su tcode.
#   Nota: los tcodes y transformaciones son IDÉNTICOS a los del script
#   mensual, porque los tcodes del dataset se definen sobre las series
#   ya en su frecuencia nativa (mensual o trimestral).
#
# ea_notransform_q()  → devuelve nivel o log-nivel (para panel E, robustez)
#   Series I(2) (tcode=5) devuelven NA → serán excluidas automáticamente.
# ============================================================

ea_transform_q <- function(x, tcode) {
  # Aplica la transformación indicada por tcode para hacer la serie estacionaria.
  # Los NAs iniciales preservan la longitud original del vector.
  if (tcode == 0)                  return(x)                          # nivel: no transformar
  if (tcode == 2 || tcode == 2.5) return(c(NA, 100 * diff(log(x))))  # tasa de crecimiento %
  if (tcode == 4 || tcode == 4.5) return(c(NA, diff(x)))             # primera diferencia
  if (tcode == 5)                  return(c(NA, NA, diff(diff(x))))   # segunda diferencia (I(2))
  warning(paste("TCODE desconocido:", tcode))
  return(x)
}

ea_notransform_q <- function(x, tcode) {
  # Devuelve la serie en su forma más natural sin transformar completamente.
  # Usado para el panel de robustez (bloque E).
  if (tcode == 0)                  return(x)               # ya en nivel
  if (tcode == 2 || tcode == 2.5) return(100 * log(x))    # log-nivel (no diferenciado)
  if (tcode == 4 || tcode == 4.5) return(x)               # nivel directo
  if (tcode == 5)                  return(rep(NA, length(x)))  # I(2) → excluir
  return(x)
}

# ============================================================
# BLOQUE C: AGREGACIÓN DE SERIES MENSUALES A TRIMESTRAL
#
# LÓGICA PASO A PASO:
#
# 1. Tomamos el dataframe mensual completo (data_raw_ea).
#
# 2. Añadimos una columna "quarter" que identifica a qué trimestre
#    pertenece cada fila:
#      Ene, Feb, Mar → 2000-01-01  (primer día del trimestre)
#      Abr, May, Jun → 2000-04-01
#      ...
#    Usamos floor_date(date, "quarter") para esto.
#
# 3. Agrupamos por trimestre y calculamos la MEDIA de los 3 meses
#    para cada serie mensual.
#
#    CONDICIÓN IMPORTANTE: na.rm = FALSE
#    Si alguno de los 3 meses tiene NA, el trimestre entero queda NA.
#    Esto evita que trimestres con datos incompletos entren al panel.
#
# 4. El resultado es un dataframe trimestral donde la fecha es el
#    primer día de cada trimestre (ej: 2000-01-01, 2000-04-01, ...).
#    Esta convención de fechas será la clave del join con las trimestrales
#    nativas en el Bloque D.
# ============================================================

cat("\n--- Agregando series mensuales a trimestral ---\n")

# Seleccionar solo las series mensuales del dataframe raw
monthly_data <- data_raw_ea %>%
  select(date, all_of(monthly_names_ea)) %>%
  mutate(quarter = floor_date(date, "quarter"))  # asignar trimestre a cada mes

# Agregar: media por trimestre para todas las series mensuales
# na.rm = FALSE → si falta algún mes, el trimestre es NA
monthly_agg <- monthly_data %>%
  group_by(quarter) %>%
  summarise(
    across(all_of(monthly_names_ea),
           ~ mean(.x, na.rm = FALSE)),
    .groups = "drop"
  ) %>%
  rename(date = quarter)

cat("Observaciones trimestrales (mensuales agregadas):", nrow(monthly_agg), "\n")
cat("Rango:", as.character(min(monthly_agg$date)),
    "a", as.character(max(monthly_agg$date)), "\n")

# ============================================================
# BLOQUE D: EXTRACCIÓN DE SERIES TRIMESTRALES NATIVAS
#
# LÓGICA PASO A PASO:
#
# 1. Las series trimestrales en el Excel tienen datos solo en el
#    ÚLTIMO MES del trimestre: marzo (mes 3), junio (6),
#    septiembre (9) y diciembre (12).
#    Los otros meses son NA.
#
# 2. Filtramos las filas del dataframe raw donde el mes es 3, 6, 9 o 12.
#    Esto nos da una fila por trimestre con los datos trimestrales.
#
# 3. Re-etiquetamos las fechas con floor_date para que coincidan con
#    la convención usada en el Bloque C (primer día del trimestre):
#      2000-03-01 → 2000-01-01  (Q1 2000)
#      2000-06-01 → 2000-04-01  (Q2 2000)
#      ...
#
# 4. Esto garantiza que ambos dataframes (mensual agregado y trimestral
#    nativo) tengan exactamente las mismas fechas para el join.
# ============================================================

cat("\n--- Extrayendo series trimestrales nativas ---\n")

quarterly_data <- data_raw_ea %>%
  select(date, all_of(quarterly_names_ea)) %>%
  filter(month(date) %in% c(3, 6, 9, 12)) %>%   # solo filas de fin de trimestre
  mutate(date = floor_date(date, "quarter"))       # re-etiquetar al inicio del trimestre

cat("Observaciones trimestrales (nativas):", nrow(quarterly_data), "\n")
cat("Rango:", as.character(min(quarterly_data$date)),
    "a", as.character(max(quarterly_data$date)), "\n")

# ============================================================
# BLOQUE E: JOIN — UNIR MENSUALES AGREGADAS + TRIMESTRALES NATIVAS
#
# Ahora ambos dataframes tienen:
#   - La misma columna "date" con fechas de inicio de trimestre
#   - Una fila por trimestre
#
# Un inner_join une solo los trimestres que existen en AMBOS datasets.
# Esto garantiza que el panel combinado no tenga trimestres con datos
# parciales de un solo lado.
#
# Resultado: data_combined_q tiene todas las series (mensuales agregadas
# + trimestrales nativas) en un único dataframe trimestral.
# ============================================================

cat("\n--- Combinando ambos panels ---\n")

data_combined_q <- inner_join(monthly_agg, quarterly_data, by = "date")

cat("Dimensiones panel combinado:", dim(data_combined_q), "\n")
cat("Rango:", as.character(min(data_combined_q$date)),
    "a", as.character(max(data_combined_q$date)), "\n")
cat("Total series disponibles:", ncol(data_combined_q) - 1, "\n")

# ============================================================
# BLOQUE F: OBJETOS PARA IVX — VERSIÓN TRIMESTRAL
#
# Construimos la regresión predictiva central en frecuencia trimestral:
#   y = inflación HICP trimestral = 100 * diff(log(HICPOV_trimestral))
#   x = tasa de desempleo UNETOT en nivel (ya está agregada en monthly_agg)
#
# DECISIÓN D2: la inflación se calcula sobre el ÍNDICE agregado (media
# de los 3 meses), no promediando inflaciones mensuales. Esto sigue la
# convención estándar del BCE y Eurostat.
#
# ALINEACIÓN TEMPORAL (igual que script mensual):
#   x[t] predice y[t+1]
#   Por eso x_full_ea_q tiene longitud T_eff + 1 (incluye el último x
#   que no tiene y correspondiente).
#   df_ea_q contiene los pares (y[t+1], x[t]) para t = 1, ..., T_eff.
# ============================================================

cat("\n--- Construyendo objetos IVX trimestrales ---\n")

# Limpiar NAs en las dos series de interés antes de diferenciar
data_ivx_q <- data_combined_q %>%
  filter(!is.na(HICPOV_EA) & !is.na(UNETOT_EA))

n_ivx_q     <- nrow(data_ivx_q)
x_full_ea_q <- data_ivx_q$UNETOT_EA                         # longitud = T_eff + 1
y_ivx_q     <- 100 * diff(log(data_ivx_q$HICPOV_EA))        # longitud = T_eff

df_ea_q <- data.frame(
  date = data_ivx_q$date[1:(n_ivx_q - 1)],
  y    = y_ivx_q,
  x    = x_full_ea_q[1:(n_ivx_q - 1)]
)

# Verificación del contrato de longitudes (condición necesaria para IVX)
stopifnot(length(x_full_ea_q) == nrow(df_ea_q) + 1)

cat("Rango df_ea_q:", as.character(min(df_ea_q$date)),
    "a", as.character(max(df_ea_q$date)), "\n")
cat("T_eff:", nrow(df_ea_q), "| x_full_ea_q:", length(x_full_ea_q), "\n")

# ============================================================
# BLOQUE G: PANEL TRANSFORMADO PARA XDLASSO — VERSIÓN TRIMESTRAL
#
# Construimos la matriz de controles W con TODAS las series del panel
# combinado (mensuales agregadas + trimestrales nativas), transformadas
# a estacionarias según sus tcodes.
#
# Variables excluidas de W (misma lógica que script mensual):
#   UNETOT_EA           → predictor de interés (d), va separado
#   UNEO25_EA, UNEU25_EA → muy correlacionadas con UNETOT_EA
#   HICPOV_EA           → variable dependiente (y)
#
# TCODES: usamos los tcodes de la columna TR1 del sheet "info".
#   Para series mensuales: tcodes_monthly
#   Para series trimestrales: tcodes_quarterly
#   La transformación se aplica DESPUÉS de la agregación — los tcodes
#   están definidos sobre series ya en frecuencia nativa.
#
# El panel resultante NO tiene desfase — main_xdlasso() hace el lag
# internamente, igual que en la versión mensual.
# ============================================================

cat("\n--- Construyendo panel XDlasso trimestral (transformado) ---\n")

vars_excluir_q <- c("UNETOT_EA", "UNEO25_EA", "UNEU25_EA", "HICPOV_EA")

# Nombres de todos los controles disponibles (mensuales + trimestrales, sin excluidas)
all_series_q <- c(monthly_names_ea, quarterly_names_ea)
vars_controles_q <- setdiff(
  intersect(all_series_q, names(data_combined_q)),  # solo las que existen en el panel
  vars_excluir_q
)

# Combinar tcodes de ambas frecuencias en un único vector
tcodes_all_q <- c(tcodes_monthly, tcodes_quarterly)

# Calcular y e x para XDlasso (sin desfasar)
infl_hicp_q     <- c(NA, 100 * diff(log(data_combined_q$HICPOV_EA)))
UNETOT_level_q  <- data_combined_q$UNETOT_EA

# Transformar cada control según su tcode
controles_trans_q <- data_combined_q %>%
  select(date, all_of(vars_controles_q))

for (v in vars_controles_q) {
  tcode_v <- tcodes_all_q[v]
  controles_trans_q[[v]] <- ea_transform_q(data_combined_q[[v]], tcode = tcode_v)
}

# Unir y y x, luego limpiar NAs
# TRNCOG_EA se excluye explícitamente porque termina antes que el resto
# del panel (igual que en 01_build_ea_data.R versión mensual).
# Sin esta exclusión, el panel se recortaría ~6 trimestres por el final.
data_trans_ea_q <- controles_trans_q %>%
  mutate(
    infl_hicp    = infl_hicp_q,
    UNETOT_level = UNETOT_level_q
  ) %>%
  select(date, UNETOT_level, infl_hicp, all_of(vars_controles_q)) %>%
  select(-any_of("TRNCOG_EA")) %>%                       # termina antes → recorta el final
  filter(!is.na(infl_hicp)) %>%                          # eliminar NA de inflación (1 obs inicial)
  filter(if_all(everything(), ~ !is.na(.)))              # eliminar cualquier columna con NAs

vars_ea_q <- setdiff(
  names(data_trans_ea_q),
  c("date", "UNETOT_level", "infl_hicp")
)

cat("Dimensiones:", dim(data_trans_ea_q), "\n")
cat("Fechas:", as.character(data_trans_ea_q$date[1]),
    "a", as.character(data_trans_ea_q$date[nrow(data_trans_ea_q)]), "\n")
cat("Controles W:", length(vars_ea_q), "\n")
cat("NAs totales:", sum(is.na(data_trans_ea_q)), "\n")

# ============================================================
# BLOQUE H: PANEL SIN TRANSFORMAR PARA XDLASSO (robustez)
#
# Igual que el Bloque E del script mensual: construimos el mismo panel
# pero con las series en nivel o log-nivel en vez de diferenciadas.
#
# Las series I(2) se excluyen directamente (tcode = 5) porque no tienen
# una representación en nivel estacionaria válida como controles W.
#
# Este panel sirve como ejercicio de robustez: si los resultados del
# XDlasso cambian mucho al usar niveles vs. diferencias, indica que
# la inferencia es sensible a la especificación de los controles.
# ============================================================

cat("\n--- Construyendo panel XDlasso trimestral (sin transformar) ---\n")

# Identificar I(2) en ambas frecuencias
vars_i2_monthly   <- names(tcodes_monthly[tcodes_monthly == 5])
vars_i2_quarterly <- names(tcodes_quarterly[tcodes_quarterly == 5])
vars_i2_q         <- c(vars_i2_monthly, vars_i2_quarterly)

vars_notrans_excluir_q   <- c(vars_excluir_q, vars_i2_q)
vars_notrans_controles_q <- setdiff(
  intersect(all_series_q, names(data_combined_q)),
  vars_notrans_excluir_q
)

controles_notrans_q <- data_combined_q %>%
  select(date, all_of(vars_notrans_controles_q))

for (v in vars_notrans_controles_q) {
  tcode_v <- tcodes_all_q[v]
  controles_notrans_q[[v]] <- ea_notransform_q(data_combined_q[[v]], tcode = tcode_v)
}

data_notrans_ea_q <- controles_notrans_q %>%
  mutate(
    infl_hicp    = infl_hicp_q,
    UNETOT_level = UNETOT_level_q
  ) %>%
  select(date, UNETOT_level, infl_hicp, all_of(vars_notrans_controles_q)) %>%
  select(-any_of("TRNCOG_EA")) %>%                       # misma exclusión que panel transformado
  filter(!is.na(infl_hicp)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

vars_notrans_ea_q <- setdiff(
  names(data_notrans_ea_q),
  c("date", "UNETOT_level", "infl_hicp")
)

cat("Dimensiones:", dim(data_notrans_ea_q), "\n")
cat("Fechas:", as.character(data_notrans_ea_q$date[1]),
    "a", as.character(data_notrans_ea_q$date[nrow(data_notrans_ea_q)]), "\n")
cat("Controles W:", length(vars_notrans_ea_q), "\n")
cat("Variables I(2) excluidas:", paste(vars_i2_q, collapse = ", "), "\n")
cat("NAs totales:", sum(is.na(data_notrans_ea_q)), "\n\n")

cat("✅ 01b_build_ea_data_quarterly.R OK\n")
cat("   Objetos IVX:     df_ea_q, x_full_ea_q\n")
cat("   Objetos XDlasso: data_trans_ea_q, vars_ea_q\n")
cat("                    data_notrans_ea_q, vars_notrans_ea_q\n")
