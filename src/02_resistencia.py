"""
02_resistencia.py
=================
Construccion de la superficie de resistencia al movimiento de fauna.
Reclasifica los rasters MapBiomas Bolivia (17 clases) en una superficie
continua de resistencia en escala relativa 1-95.

Parametrizacion basada en:
    - Zeller et al. (2012). Landscape Ecology 27(6): 777-797.
    - Wade et al. (2015). USDA-FS GTR RMRS-333.
    - Rabinowitz & Zeller (2010). Biol. Conserv. 143(4): 939-945.
    - Naiman et al. (2005). Riparia. Elsevier Academic Press.
    - Noss et al. (2003). Conserv. Biol. 17(2): 601-611.

Especies paraguas: jaguar (Panthera onca) y tapir de tierras bajas
    (Tapirus terrestris).

Uso:
    python 02_resistencia.py
    # Los rasters de salida se guardan en data/rasters_90m/resistencia_{anio}.tif

Requisitos:
    numpy>=1.26, gdal>=3.12
    Entorno: conda activate conectividad-pirai

Autor: Gloria Eliana Torrez Castro — PPGG/UFC, 2024-2025
"""

import numpy as np
import os
from osgeo import gdal, gdalconst

# ── Rutas ───────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_DIR   = os.path.join(BASE_DIR, "data", "rasters_90m")   # rasters MapBiomas 90m EPSG:32720
OUT_DIR  = os.path.join(BASE_DIR, "data", "rasters_90m")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Tabla de resistencia — 17 clases MapBiomas Bolivia Col. 2024 ─────────────
# Formato: codigo_mapbiomas -> valor_resistencia (1-95)
# Escala: 1 = maxima permeabilidad (habitat fuente)
#         95 = maxima resistencia (barrera total)
#
# Fuente parametrizacion: ver docstring del modulo.
# Archivo CSV de referencia: data/resistance_table.csv
TABLA_RESISTENCIA = {
    # ── Habitat fuente ─────────────────────────────────────────────────────
    3:  1,    # Formacion bocosa — bosque nativo maduro
    6:  1,    # Bosque inundable
    68: 1,    # Bosque nativo secundario
    72: 1,    # Formacion bocosa en general
    # ── Corredor ripario ───────────────────────────────────────────────────
    11: 3,    # Humedal / campo humedo (eje de conectividad fluvial)
    33: 3,    # Cuerpo de agua / rio principal
    # ── Matriz semipermeable ───────────────────────────────────────────────
    12: 15,   # Pastizal natural (campo limpio)
    13: 20,   # Formacion herbacea y arbustiva
    4:  22,   # Formacion savana cerrado (Chiquitania)
    # ── Barrera moderada ───────────────────────────────────────────────────
    29: 35,   # Afloramiento rocoso
    31: 40,   # Acuicultura / estanques artificiales
    9:  50,   # Silvicultura (plantaciones forestales)
    # ── Barrera alta ───────────────────────────────────────────────────────
    15: 55,   # Pasturas manejadas (riesgo de caza)
    21: 60,   # Mosaico agropecuario
    18: 65,   # Agricultura — cultivos anuales y permanentes
    48: 65,   # Otros cultivos no identificados
    36: 70,   # Arrozal / cultivo irrigado (inundacion periodica)
    # ── Barrera severa ─────────────────────────────────────────────────────
    25: 85,   # Suelo desnudo / area quemada
    24: 90,   # Area urbanizada / infraestructura
    # ── Barrera maxima ─────────────────────────────────────────────────────
    30: 95,   # Mineria (disturbio permanente)
}

# Valor por defecto para codigos no listados (clases menores o no cartografiadas)
RESISTENCIA_DEFAULT = 40

# Valor de resistencia para NoData (se mantiene como NoData en la salida)
NODATA_OUT = -9999


def construir_resistencia(anio: int) -> str:
    """
    Reclasifica el raster MapBiomas de un anio en superficie de resistencia.

    Args:
        anio: Anio de analisis.

    Returns:
        Ruta del raster de resistencia generado.
    """
    ruta_in  = os.path.join(IN_DIR, f"mb_riogrande_{anio}_90m.tif")
    ruta_out = os.path.join(OUT_DIR, f"resistencia_{anio}_90m.tif")

    ds = gdal.Open(ruta_in, gdalconst.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f"No se encontro: {ruta_in}")

    band = ds.GetRasterBand(1)
    arr  = band.ReadAsArray().astype(np.int16)
    nd   = int(band.GetNoDataValue() or -9999)
    gt   = ds.GetGeoTransform()
    proj = ds.GetProjection()
    rows, cols = arr.shape
    ds = None

    # Reclasificar
    res_arr = np.full_like(arr, RESISTENCIA_DEFAULT, dtype=np.int16)
    for cod, val in TABLA_RESISTENCIA.items():
        res_arr[arr == cod] = val
    res_arr[arr == nd] = NODATA_OUT   # preservar NoData

    # Escribir raster de salida
    driver  = gdal.GetDriverByName("GTiff")
    ds_out  = driver.Create(ruta_out, cols, rows, 1, gdal.GDT_Int16,
                            ["COMPRESS=LZW", "TILED=YES"])
    ds_out.SetGeoTransform(gt)
    ds_out.SetProjection(proj)
    b_out = ds_out.GetRasterBand(1)
    b_out.SetNoDataValue(NODATA_OUT)
    b_out.WriteArray(res_arr)
    ds_out.FlushCache()
    ds_out = None

    print(f"  {anio}: resistencia guardada -> {os.path.basename(ruta_out)}")
    return ruta_out


# ── Ejecutar ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Para Circuitscape se necesitan los 3 cortes de grafos
    ANIOS_CS = [1985, 2005, 2024]

    print("=" * 65)
    print("SUPERFICIE DE RESISTENCIA — Subcuenca rio Grande-Pirai")
    print("=" * 65)
    for anio in ANIOS_CS:
        construir_resistencia(anio)
    print("\nOK: superficies de resistencia generadas para Circuitscape")
    print("    Siguiente paso: ejecutar 03_nodos.py para generar los")
    print("    rasters de nodos nucleo, o ejecutar Circuitscape con")
    print("    los archivos .ini de /circuitscape/")
