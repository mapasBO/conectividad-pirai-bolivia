# Conectividad ecológica — Subcuenca río Grande–Piraí, Bolivia (1985–2024)

**Artículo:** *Colapso silencioso: pérdida no lineal de conectividad funcional (−35,6%) frente a pérdida moderada de área forestal (−13,2%) en la subcuenca del río Grande–Piraí, Bolivia, 1985–2024*

**Autora:** Glória Eliana Torrez Castro — PPGG/UFC, 2024–2025  
**Zenodo:** https://doi.org/10.5281/zenodo.XXXXXXX *(DOI a asignar en revisión final)*

---

## Descripción

Código fuente, archivos de configuración y tabla de parámetros para el análisis
multitemporal de conectividad funcional del paisaje forestal en la subcuenca del
río Grande–Piraí (Santa Cruz, Bolivia), período 1985–2024.

El flujo reproduce el proceso estandarizado de modelado de conectividad descrito
por Wade et al. (2015) e implementa los índices Conefor Sensinode (Saura & Torné,
2009) mediante NetworkX en Python.

---

## Estructura del repositorio

```
conectividad-pirai-bolivia/
├── src/
│   ├── 01_fragmentacion.py   # Métricas FRAGSTATS (NP, CA, LPI, MPS, FRAC, ENN, SHDI)
│   ├── 02_resistencia.py     # Superficie de resistencia (17 clases, escala 1-95)
│   ├── 03_nodos.py           # Extracción de nodos núcleo ≥ 500 ha → ASC para Circuitscape
│   ├── 04_conectividad.py    # Índices IIC, PC y dPC con NetworkX
│   └── 05_proyeccion.py      # Proyección PC → 2040 bajo 3 escenarios
├── circuitscape/
│   ├── config_1985_90m.ini   # Configuración Circuitscape año 1985
│   ├── config_2005_90m.ini   # Configuración Circuitscape año 2005
│   └── config_2024_90m.ini   # Configuración Circuitscape año 2024
├── data/
│   ├── resistance_table.csv  # Tabla de resistencia: 17 clases + fuentes bibliográficas
│   ├── rasters_30m/          # Rásteres MapBiomas originales 30 m (Zenodo)
│   ├── rasters_90m/          # Rásteres remuestreados 90 m EPSG:32720 (Zenodo)
│   ├── fragmentacion/        # CSVs con métricas FRAGSTATS por año
│   └── conectividad/         # CSVs con índices PC, IIC, dPC y proyecciones
├── environment.yml           # Entorno conda reproducible
└── README.md
```

---

## Instalación rápida

```bash
# 1. Clonar el repositorio
git clone https://github.com/mapasBO/conectividad-pirai-bolivia.git
cd conectividad-pirai-bolivia

# 2. Crear entorno conda
conda env create -f environment.yml
conda activate conectividad-pirai

# 3. Descargar datos desde Zenodo (rásteres preprocesados)
#    https://doi.org/10.5281/zenodo.XXXXXXX
#    Colocar en data/rasters_30m/ y data/rasters_90m/
```

---

## Flujo de ejecución

```bash
# Paso 1: Métricas de fragmentación (7 años × 5 clases)
python src/01_fragmentacion.py
# Salida: data/fragmentacion/fragmentacion_riogrande_1985_2024.csv
# Tiempo: ~15 min

# Paso 2: Superficie de resistencia (3 años para Circuitscape)
python src/02_resistencia.py
# Salida: data/rasters_90m/resistencia_{1985,2005,2024}_90m.tif
# Tiempo: ~3 min

# Paso 3: Extracción de nodos núcleo
python src/03_nodos.py
# Salida: data/rasters_90m/nodos_{anio}_90m.asc
#         data/fragmentacion/nodos_atributos_{anio}.csv
# Tiempo: ~5 min

# Paso 4a: Ejecutar Circuitscape (GUI o línea de comandos)
#   Windows: Circuitscape.exe config_1985_90m.ini
#   Python:  circuitscape_compute('circuitscape/config_1985_90m.ini')
# Repetir para 2005 y 2024. Tiempo: 5-15 min por año.

# Paso 4b: Índices de conectividad (IIC, PC, dPC)
python src/04_conectividad.py
# Salida: data/conectividad/indices_PC_IIC_dPC.csv
# Tiempo: ~30-60 min (48 iteraciones dPC)

# Paso 5: Proyección PC → 2040
python src/05_proyeccion.py
# Salida: data/conectividad/proyeccion_PC_2040_escenarios.csv
# Tiempo: < 1 min
```

---

## Parámetros clave

| Parámetro | Valor | Fundamento |
|-----------|-------|------------|
| Resolución Circuitscape | 90 m | Balance cómputo/precisión |
| Vecindad | 8 celdas | Conectividad diagonal activa |
| Umbral nodo | ≥ 500 ha | Rabinowitz & Zeller (2010) |
| Umbral conexión | 5 km | Distancia Euclidiana entre centroides |
| Distancia dispersión efectiva | 3 km | Noss et al. (2003); Zeller et al. (2012) |
| Función dispersión | p(d) = exp(−d/3) | Decaimiento exponencial negativo |
| Área del paisaje | 1.352.700 ha | Subcuenca río Grande–Piraí |
| Umbral inviabilidad | PC < 0,020 | Bennett (2003) |
| Especies paraguas | Jaguar + tapir | Rabinowitz & Zeller (2010); Noss et al. (2003) |

---

## Requisitos de hardware

- RAM: 16 GB mínimo (recomendado 32 GB para Circuitscape con 16 nodos a 90 m)
- CPU: 4 núcleos mínimo
- Almacenamiento: ~8 GB para datos y resultados
- Tiempo total estimado: 4–6 horas

---

## Herramientas y versiones

| Herramienta | Versión | Licencia |
|-------------|---------|----------|
| Python | 3.12 | PSF |
| QGIS | 3.44.8-Solothurn | GPL-2.0 |
| GDAL/OGR | 3.12.2 | MIT/X |
| NumPy | 1.26 | BSD-3 |
| SciPy | 1.11 | BSD-3 |
| NetworkX | 3.2.1 | BSD-3 |
| Circuitscape | 4.0.5 | MIT |
| MapBiomas Bolivia | Colección 2024 | CC BY 4.0 |

---

## Cita

```
Torrez Castro, G.E. (2025). Colapso silencioso: pérdida no lineal de
conectividad funcional en la subcuenca del río Grande–Piraí, Bolivia,
1985–2024. PPGG/UFC. https://doi.org/10.5281/zenodo.XXXXXXX
```

---

## Referencias principales

- Bennett, A.F. (2003). *Linkages in the landscape* (2ª ed.). IUCN.
- McRae et al. (2008). Ecology 89(10): 2712–2724. https://doi.org/10.1890/07-1861.1
- Saura & Torné (2009). Env. Modelling & Software 24(1): 135–139.
- Wade et al. (2015). USDA-FS GTR RMRS-333.
- Zeller et al. (2012). Landscape Ecology 27(6): 777–797.
