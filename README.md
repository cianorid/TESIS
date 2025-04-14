# 📊 Tesis: Modelos de Volatilidad Estocástica y Filtro de Partículas

Este repositorio contiene el código desarrollado como parte de una tesis de grado centrada en el modelado de la volatilidad en series temporales financieras. El objetivo principal es implementar un modelo de **volatilidad estocástica (SV)** y aplicar técnicas de **inferencia bayesiana (MCMC)** y **filtros secuenciales (filtro de partículas)** para estimar la volatilidad latente a partir de datos reales del mercado.

---

## 🧠 Contenido

El proyecto está dividido en los siguientes scripts:

| Archivo                      | Descripción                                                                 |
|-----------------------------|-----------------------------------------------------------------------------|
| `TESIS_simulacion.R`        | Simulación de datos bajo un modelo SV y evaluación con parámetros conocidos. |
| `TESIS_funciones.R`         | Funciones auxiliares para el algoritmo MCMC: actualizaciones de parámetros. |
| `TESIS_filtro_de_particulas.R` | Implementación del filtro de partículas para estimar la volatilidad.          |
| `TESIS_datos.R`             | Análisis completo con datos reales de Bitcoin, incluyendo diagnóstico y comparación empírica de distribuciones. |

---

## 🔧 Requisitos

Para correr los scripts se necesitan las siguientes librerías de R:

- `MASS`
- `coda`
- `ggplot2`
- `quantmod`
- `gridExtra`
- `plotly`

Podés instalarlas con:

```r
install.packages(c("MASS", "coda", "ggplot2", "quantmod", "gridExtra", "plotly"))
