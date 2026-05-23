install.packages("readxl")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("vars")
install.packages("here")
install.packages("tidyr")
install.packages("lubridate")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("gridExtra")
install.packages("gt")
library(MASS)
library(gt)
library(lubridate)
library(here)
library(tidyr)
library(vars)
library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)
library(gridExtra)



# EXPORTO DATA
GDP   <- read_excel(here("data", "GDP - USA.xlsx"))
HOURS <- read_excel(here("data", "HOURS - USA.xlsx"))
empleo_mensual   <- read_excel(here("data", "EMP - MONTHLY - USA.xlsx"))

# EMPLEO TRIMESTRAL 

EMPLOYMENT <- empleo_mensual %>%
  # Asegurar que 'date' sea reconocido como fecha
  mutate(observation_date = as.Date(observation_date)) %>% 
  # Agrupar por el inicio de cada trimestre (enero, abril, julio, octubre)
  group_by(observation_date = floor_date(observation_date, "quarter")) %>%
  # Calcular el promedio de los 3 meses del trimestre
  summarise(CE16OV = mean(CE16OV, na.rm = TRUE)) %>%
  ungroup()


# RENOMBRO LAS COLUMNAS
GDP <- GDP %>%
  rename(GDPR = GDPC1)

EMPLOYMENT <- EMPLOYMENT %>%
  rename(EMP = CE16OV)

HOURS <- HOURS %>%
  rename(HRS = HOANBS)

# DELIMITO 1948:1 A 1994:4

GDP_cleaned <- GDP %>%
  mutate(observation_date = as.Date(observation_date)) %>% 
  filter(observation_date >= as.Date("1948-01-01") & 
           observation_date <= as.Date("1994-10-01"))

EMP_cleaned <- EMPLOYMENT %>%
  mutate(observation_date = as.Date(observation_date)) %>% 
  filter(observation_date >= as.Date("1948-01-01") & 
           observation_date <= as.Date("1994-10-01"))

HRS_cleaned <- HOURS %>%
  mutate(observation_date = as.Date(observation_date)) %>% 
  filter(observation_date >= as.Date("1948-01-01") & 
           observation_date <= as.Date("1994-10-01"))

#APLICAMOS LOGARITMO 

GDP_cleaned <- GDP_cleaned %>%
  mutate(y = log(GDPR))

EMP_cleaned <- EMP_cleaned %>%
  mutate(e = log(EMP))

HRS_cleaned <- HRS_cleaned %>%
  mutate(n = log(HRS))

# PRODUCTIVIDAD HRS y PRODUCTIVIDAD EMP

data_final <- GDP_cleaned %>%
  inner_join(HRS_cleaned, by = "observation_date") %>%
  inner_join(EMP_cleaned, by = "observation_date") 

productivity_table <- data_final %>%
  mutate(x_hrs = y - n) %>%
  mutate(x_emp = y - e) %>%
  select(observation_date, y, n, e, x_hrs, x_emp)

# BASE DE DATOS FINALES 
DatosFinal <- productivity_table %>%
  rename(date = observation_date) 

# Tasas de crecimiento 

DatosFinal <- DatosFinal %>%
  # Aseguramos el orden cronológico
  arrange(date) %>% 
  mutate(
    # Tasa de crecimiento de la productividad (horas y empleo)
    dx_hrs = x_hrs - lag(x_hrs),
    dx_emp = x_emp - lag(x_emp),
    # Tasa de crecimiento del insumo de trabajo (horas y empleo)
    dn = n - lag(n),
    de = e - lag(e)
  )

#######
# TABLE 1 
# PANEL A

# Filtramos NA causados por la primera diferencia
df_hrs <- na.omit(DatosFinal[, c("dx_hrs", "dn")])
df_emp <- na.omit(DatosFinal[, c("dx_emp", "de")])


# FUNCIÓN AUXILIAR PARA CALCULAR CORRELACIONES (Tabla 1 - Panel A)

calc_correlations <- function(df, var_names) {
  
  # A. Correlación Incondicional (Directamente de los datos)
  unc_cor <- cor(df[,1], df[,2])
  
  # B. Estimación del VAR Reducido (Usando 4 rezagos)
  # Se utiliza "const" para incluir constante en la regresión
  var_est <- VAR(df, p = 4, type = "const")
  
  # C. Estimación del SVAR con Restricción de Blanchard-Quah (Largo Plazo) (EC. 23)
  # La función BQ requiere que la variable que asume el impacto permanente esté primero.
  # MATRIZ LOWER TRIANGULAR SEGUN RETRICCION.
  svar_bq <- BQ(var_est)
  
  # D. Extracción de Coeficientes de Impulso-Respuesta Estructurales (MA)
  # Limitamos la suma infinita a j = 100 periodos (suficiente para convergencia)
  # Obtiene C(L)
  irfs <- irf(svar_bq, n.ahead = 100, boot = FALSE)
  
  # C_tech = Respuestas al shock 1 (Tecnológico)
  # C_nontech = Respuestas al shock 2 (No Tecnológico)
  C_tech <- irfs$irf[[1]] 
  C_nontech <- irfs$irf[[2]]
  
  # E. Cálculo de Correlación Condicional al Shock Tecnológico (Ec. 24)
  cov_z <- sum(C_tech[, var_names[1]] * C_tech[, var_names[2]])
  var_x_z <- sum(C_tech[, var_names[1]]^2)
  var_n_z <- sum(C_tech[, var_names[2]]^2)
  cond_cor_tech <- cov_z / sqrt(var_x_z * var_n_z) # Ecuación 24 para choque tec
  
  # F. Cálculo de Correlación Condicional al Shock No Tecnológico (Ec. 24)
  cov_m <- sum(C_nontech[, var_names[1]] * C_nontech[, var_names[2]])
  var_x_m <- sum(C_nontech[, var_names[1]]^2)
  var_n_m <- sum(C_nontech[, var_names[2]]^2)
  cond_cor_nontech <- cov_m / sqrt(var_x_m * var_n_m) # Ecuación 24 para choque no tec
  
  # Retornar vector de resultados
  return(c(Unconditional = unc_cor, 
           Tech_Conditional = cond_cor_tech, 
           NonTech_Conditional = cond_cor_nontech))
}


# EJECUCIÓN DEL MODELO PARA HORAS Y EMPLEO

# Resultados para Horas (Panel A - Fila 1)
res_hours <- calc_correlations(df_hrs, c("dx_hrs", "dn"))

# Resultados para Empleo (Panel A - Fila 2)
res_employment <- calc_correlations(df_emp, c("dx_emp", "de"))

# Compilar la Tabla
Tabla1_PanelA <- rbind(Hours = res_hours, Employment = res_employment)
print(round(Tabla1_PanelA, 2))

###### FIGURA 1: PRODUCTIVDAD vs HORAS: DATA & CHOQUES

# multiplico por 100 los datos ## reescalado
df_hrs <- df_hrs*100
#########

# 0. Definimos en un ambiente general 
var_est <- VAR(df_hrs, p = 4, type = "const") # se estima el modelo VAR
svar_bq <- BQ(var_est) # se aplica la restricción estructural

# 1. Recuperar los componentes históricos
# Extraemos los residuos del VAR estimado previamente (var_res)
residuos <- residuals(var_est) # extrae residuos (ut) del modelo etimado
S_mat <- summary(svar_bq)$A # Extrae la matriz de impacto contemporaneo
# esta matriz mapea como los choques estructurales afectan a las variables en t
shocks_e <- t(solve(S_mat) %*% t(residuos)) # recupera choques estructurales (Et) ortogonales. 

# 2. Reconstruir los componentes (Descomposición Histórica)
# Proyectamos los shocks sobre la estructura del modelo 
n_obs <- nrow(residuos)
tech_data <- matrix(0, n_obs, 2)
non_tech_data <- matrix(0, n_obs, 2)

# Usamos la representación de Media Móvil (MA) 
Phi <- irf(svar_bq, n.ahead = n_obs)$irf # reccupera los coeficientes de MEDIAS MOVILES
# Estos coeficientes Cj dictan la respuesta de las variables a cada choque a lo largo del tiempo

#esta es la descomposicion historica de Wold. 
for(i in 1:n_obs) {
  for(j in 0:(i-1)) {
    # Componente Tecnológico (Shock 1) [cite: 274, 275]
    tech_data[i,] <- tech_data[i,] + Phi[[1]][j+1,] * shocks_e[i-j, 1] #componente impulsado por choques tec
    # Componente No Tecnológico (Shock 2) [cite: 274, 275]
    non_tech_data[i,] <- non_tech_data[i,] + Phi[[2]][j+1,] * shocks_e[i-j, 2] # choques no tec
  }
}

# aquí se grafican los residuos, pero yo debo graficar los las variables originales.
# 3. Creación de las 3 nubes de puntos (Figure 1 del paper) 
df_plot <- data.frame(
  dn = residuos[,2], dx = residuos[,1],
  tech_dn = tech_data[,2], tech_dx = tech_data[,1],
  non_dn = non_tech_data[,2], non_dx = non_tech_data[,1]
)

############### reescalada
# Definimos un estilo de cuadrícula común para las tres gráficas
estilo_gali <- theme_light() + 
  theme(
    panel.grid.major = element_line(color = "gray90"), # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray95"), # Cuadrícula secundaria
    panel.background = element_rect(fill = "white", color = "black"), # Fondo blanco y borde negro
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14), # Título centrado
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )

# Gráfica A: Data
g1 <- ggplot(df_plot, aes(x=dn, y=dx)) + 
  geom_point(shape=1, alpha=0.7) + 
  geom_smooth(method="lm", color="black", se=FALSE, linewidth=0.5) + 
  labs(title="Data", x="hours", y="productivity") + 
  scale_x_continuous(limits = c(-3, 4)) +
  scale_y_continuous(limits = c(-1.5, 3.0)) +
  estilo_gali

# Gráfica B: Technology Component
g2 <- ggplot(df_plot, aes(x=tech_dn, y=tech_dx)) + 
  geom_point(shape=1, alpha=0.7) + 
  geom_smooth(method="lm", color="black", se=FALSE, linewidth=0.5) + 
  labs(title="Technology Component", x="hours", y="productivity") + 
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2.0, 2.5)) +
  estilo_gali

# Gráfica C: Nontechnology Component
g3 <- ggplot(df_plot, aes(x=non_dn, y=non_dx)) + 
  geom_point(shape=1, alpha=0.7) + 
  geom_smooth(method="lm", color="black", se=FALSE, linewidth=0.5) + 
  labs(title="Nontechnology Component", x="hours", y="productivity") + 
  scale_x_continuous(limits = c(-2.8, 2.1)) +
  scale_y_continuous(limits = c(-1.05, 1.40)) +
  estilo_gali

# Mostrar las gráficas combinadas
grid.arrange(g1, g2, g3, ncol=1)



############# otro diseño
# Gráfica A: Data (Correlación Incondicional ~ -0.26) [cite: 305, 325]
g1 <- ggplot(df_plot, aes(x=dn, y=dx)) + geom_point(alpha=0.5) + 
  geom_smooth(method="lm", color="black", se=F) + labs(title="Data", x="hours", y="productivity") + theme_minimal()

# Gráfica B: Technology Component (Correlación Condicional ~ -0.82) [cite: 310, 325]
g2 <- ggplot(df_plot, aes(x=tech_dn, y=tech_dx)) + geom_point(alpha=0.5) + 
  geom_smooth(method="lm", color="red", se=F) + labs(title="Technology Component", x="hours", y="productivity") + theme_minimal()

# Gráfica C: Nontechnology Component (Correlación Condicional ~ 0.26) [cite: 310, 325, 341]
g3 <- ggplot(df_plot, aes(x=non_dn, y=non_dx)) + geom_point(alpha=0.5) + 
  geom_smooth(method="lm", color="blue", se=F) + labs(title="Nontechnology Component", x="hours", y="productivity") + theme_minimal()

# Mostrar resultados
grid.arrange(g1, g2, g3, ncol=1)

#### FIGURA 2: IMPULSO RESPUESTA

# 2. Re-estimar el VAR bivariado base (Panel A: primeras diferencias)
var_base <- VAR(df_hrs, p = 4, type = "const")
svar_bq_base <- BQ(var_base)

# 3. Calcular Funciones Impulso-Respuesta ACUMULADAS (para ver el efecto en niveles)
# n.ahead = 12 (Galí grafica hasta 12 trimestres)
# cumulative = TRUE (Suma los impactos para pasar de diferencia a nivel)
# boot = TRUE (Calcula los intervalos de confianza, tardará unos segundos)
irf_acum <- irf(svar_bq_base, n.ahead = 12, cumulative = TRUE, boot = TRUE, ci = 0.95)

# 4. Extraer los datos para las gráficas
# Shock 1 = Tecnología | Shock 2 = No Tecnología
# Variable 1 = Productividad (dx_hrs) | Variable 2 = Horas (dn)

# --- SHOCK TECNOLÓGICO ---
tech_prod <- irf_acum$irf[[1]][, 1]
tech_prod_low <- irf_acum$Lower[[1]][, 1]
tech_prod_up <- irf_acum$Upper[[1]][, 1]

tech_hours <- irf_acum$irf[[1]][, 2]
tech_hours_low <- irf_acum$Lower[[1]][, 2]
tech_hours_up <- irf_acum$Upper[[1]][, 2]

# Recuperar el PIB (y = x + n)
tech_gdp <- tech_prod + tech_hours
tech_gdp_low <- tech_prod_low + tech_hours_low # Aproximación lineal de bandas
tech_gdp_up <- tech_prod_up + tech_hours_up

# --- SHOCK NO TECNOLÓGICO (Demanda) ---
nontech_prod <- irf_acum$irf[[2]][, 1]
nontech_prod_low <- irf_acum$Lower[[2]][, 1]
nontech_prod_up <- irf_acum$Upper[[2]][, 1]

nontech_hours <- irf_acum$irf[[2]][, 2]
nontech_hours_low <- irf_acum$Lower[[2]][, 2]
nontech_hours_up <- irf_acum$Upper[[2]][, 2]

# Recuperar el PIB (y = x + n)
nontech_gdp <- nontech_prod + nontech_hours
nontech_gdp_low <- nontech_prod_low + nontech_hours_low
nontech_gdp_up <- nontech_prod_up + nontech_hours_up

# 5. Configurar el lienzo para 6 gráficas (3 filas x 2 columnas)
par(mfrow = c(3, 2), mar = c(3, 4, 3, 1), oma = c(2, 0, 2, 0))

# Función auxiliar para dibujar cada gráfica con el estilo de Galí
plot_gali <- function(pe, lower, upper, title, ylab) {
  plot(0:12, pe, type = "l", lwd = 2, ylim = range(c(lower, upper, 0, pe)),
       xlab = "", ylab = ylab, main = title, cex.main = 1.2)
  lines(0:12, lower, lty = 2, col = "red") # Banda inferior
  lines(0:12, upper, lty = 2, col = "red") # Banda superior
  abline(h = 0, col = "black", lwd = 1)    # Línea cero
  points(0:12, pe, pch = 17, cex = 1)      # Triángulos como en el paper original
}

# --- Fila 1: Productividad ---
plot_gali(tech_prod, tech_prod_low, tech_prod_up, "Technology Shock", "productivity")
plot_gali(nontech_prod, nontech_prod_low, nontech_prod_up, "Nontechnology Shock", "")

# --- Fila 2: PIB (GDP) ---
plot_gali(tech_gdp, tech_gdp_low, tech_gdp_up, "", "gdp")
plot_gali(nontech_gdp, nontech_gdp_low, nontech_gdp_up, "", "")

# --- Fila 3: Horas ---
plot_gali(tech_hours, tech_hours_low, tech_hours_up, "", "hours")
plot_gali(nontech_hours, nontech_hours_low, nontech_hours_up, "", "")

# Título
mtext("Trimestres tras el shock", side = 1, outer = TRUE, cex = 0.9)

##### prueba para errores estandar ######## 

calc_correlations_with_significance <- function(df, var_names, lags = 4, horizon = 100, draws = 500) {
  
  # -------------------------------------------------------------------
  # 1. ESTIMACIÓN PUNTUAL (MODELO ORIGINAL)
  # -------------------------------------------------------------------
  var_original <- VAR(df, p = lags, type = "const")
  svar_bq <- BQ(var_original)
  irfs_orig <- irf(svar_bq, n.ahead = horizon, boot = FALSE)
  
  compute_conditional_corrs <- function(irf_object) {
    # Componente Tecnológico
    C_tech <- irf_object$irf[[1]]
    cov_z <- sum(C_tech[, var_names[1]] * C_tech[, var_names[2]])
    var_x_z <- sum(C_tech[, var_names[1]]^2)
    var_n_z <- sum(C_tech[, var_names[2]]^2)
    cond_cor_tech <- cov_z / sqrt(var_x_z * var_n_z)
    
    # Componente No Tecnológico
    num_shocks <- length(irf_object$irf)
    var_x_m <- 0
    var_n_m <- 0
    cov_m <- 0
    
    for (i in 2:num_shocks) {
      C_nontech <- irf_object$irf[[i]]
      var_x_m <- var_x_m + sum(C_nontech[, var_names[1]]^2)
      var_n_m <- var_n_m + sum(C_nontech[, var_names[2]]^2)
      cov_m <- cov_m + sum(C_nontech[, var_names[1]] * C_nontech[, var_names[2]])
    }
    cond_cor_nontech <- cov_m / sqrt(var_x_m * var_n_m)
    
    return(c(Tech = cond_cor_tech, NonTech = cond_cor_nontech))
  }
  
  point_estimates <- compute_conditional_corrs(irfs_orig)
  
  # -------------------------------------------------------------------
  # 2. PROCESO DE MONTE CARLO PARA ERRORES ESTÁNDAR
  # -------------------------------------------------------------------
  B_mat <- Bcoef(var_original)
  Sigma_u <- summary(var_original)$covres
  
  K <- ncol(df)
  T_len <- nrow(df)
  init_Y <- as.matrix(df[1:lags, ])
  
  cor_tech_draws <- numeric(draws)
  cor_nontech_draws <- numeric(draws)
  
  set.seed(1999) # Replicabilidad
  valid_draws <- 0
  max_attempts <- draws * 5
  attempts <- 0
  
  while (valid_draws < draws) {
    attempts <- attempts + 1
    u_sim <- MASS::mvrnorm(n = T_len, mu = rep(0, K), Sigma = Sigma_u)
    
    Y_sim <- matrix(0, nrow = T_len, ncol = K)
    colnames(Y_sim) <- colnames(df)
    Y_sim[1:lags, ] <- init_Y
    
    const_vec <- B_mat[, ncol(B_mat)]
    A_matrices <- B_mat[, 1:(K * lags)]
    
    for (t in (lags + 1):T_len) {
      y_t <- const_vec + u_sim[t, ]
      for (j in 1:lags) {
        A_j <- A_matrices[, ((j - 1) * K + 1):(j * K)]
        y_t <- y_t + A_j %*% Y_sim[t - j, ]
      }
      Y_sim[t, ] <- y_t
    }
    
    draw_result <- tryCatch({
      var_sim <- VAR(as.data.frame(Y_sim), p = lags, type = "const")
      compute_conditional_corrs(irf(BQ(var_sim), n.ahead = horizon, boot = FALSE))
    }, error = function(e) NULL)
    
    if (!is.null(draw_result) && !any(is.na(draw_result))) {
      valid_draws <- valid_draws + 1
      cor_tech_draws[valid_draws] <- draw_result["Tech"]
      cor_nontech_draws[valid_draws] <- draw_result["NonTech"]
    }
    
    if (attempts > max_attempts) stop("Inestabilidad excesiva en el Monte Carlo.")
  }
  
  se_tech <- sd(cor_tech_draws)
  se_nontech <- sd(cor_nontech_draws)
  
  # -------------------------------------------------------------------
  # 3. CÁLCULO DE SIGNIFICANCIA ESTADÍSTICA (PRUEBA Z DE 2 COLAS)
  # -------------------------------------------------------------------
  # Calculamos el estadístico z asintótico
  z_tech <- point_estimates["Tech"] / se_tech
  z_nontech <- point_estimates["NonTech"] / se_nontech
  
  # Calculamos p-valores de dos colas
  pval_tech <- 2 * (1 - pnorm(abs(z_tech)))
  pval_nontech <- 2 * (1 - pnorm(abs(z_nontech)))
  
  # Asignamos asteriscos según los umbrales de Galí
  get_asterisks <- function(pval) {
    if (pval <= 0.05) return("**")
    if (pval <= 0.10) return("*")
    return("")
  }
  
  ast_tech <- get_asterisks(pval_tech)
  ast_nontech <- get_asterisks(pval_nontech)
  
  # Formateamos la salida para que coincida visualmente con el paper
  formatted_tech <- sprintf("%.2f%s (%.2f)", point_estimates["Tech"], ast_tech, se_tech)
  formatted_nontech <- sprintf("%.2f%s (%.2f)", point_estimates["NonTech"], ast_nontech, se_nontech)
  
  res_table <- data.frame(
    Estimate = c(point_estimates["Tech"], point_estimates["NonTech"]),
    Std_Error = c(se_tech, se_nontech),
    Z_stat = c(z_tech, z_nontech),
    P_value = c(pval_tech, pval_nontech),
    Significance = c(ast_tech, ast_nontech),
    Gali_Format = c(formatted_tech, formatted_nontech),
    row.names = c("Technology_Conditional", "Nontechnology_Conditional")
  )
  
  return(res_table)
}

# =====================================================================
# EJECUCIÓN Y PRESENTACIÓN DE RESULTADOS
# =====================================================================
cat("=== Panel A: Growth rates (Hours - Multivariate) ===\n")
res_hours <- calc_correlations_with_significance(df_hrs, var_names = c("dx_hrs", "dn"), draws = 500)

# Mostramos la tabla completa con métricas detalladas
print(res_hours[, c("Estimate", "Std_Error", "P_value", "Gali_Format")])
cat("\n")

cat("=== Panel A: Growth rates (Employment - Multivariate) ===\n")
res_employment <- calc_correlations_with_significance(df_emp, var_names = c("dx_emp", "de"), draws = 500)
print(res_employment[, c("Estimate", "Std_Error", "P_value", "Gali_Format")])


