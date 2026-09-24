"""
03_nodos.py
===========
Extraccion de nodos nucleo de habitat para Circuitscape y NetworkX.
Identifica parches de bosque >= 500 ha y genera:
    - Raster de nodos para Circuitscape (.asc, un ID por parche)
    - CSV de atributos de nodos (ID, area, centroide, bbox)
    - Shapefile de nodos (opcional, requiere QGIS o geopandas)

Umbral de nodo: >= 500 ha, fundamentado en los requerimientos minimos
de territorio para jaguares adultos en el Chaco boliviano
(Rabinowitz & Zeller, 2010) y en la sensibilidad del tapir a la
fragmentacion en ecosistemas de transicion Chaco-Chiquitania
(Noss et al., 2003).

Uso:
    python 03_nodos.py

Requisitos:
    numpy>=1.26, scipy>=1.11, gdal>=3.12
    Entorno: conda activate conectividad-pirai

Autor: Gloria Eliana Torrez Castro — PPGG/UFC, 2024-2025
"""

import numpy as np
import os
import csv
from scipy import ndimage
from osgeo import gdal, osr

# ── Rutas ───────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_DIR   = os.path.join(BASE_DIR, "data", "rasters_90m")
OUT_DIR  = os.path.join(BASE_DIR, "data", "rasters_90m")
CSV_DIR  = os.path.join(BASE_DIR, "data", "fragmentacion")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parametros ───────────────────────────────────────────────────────────────
# Codigos MapBiomas Bolivia Col. 2024 que corresponden a Bosque (hábitat fuente)
CODIGOS_BOSQUE = [3, 6, 11, 68, 72]

# Umbral minimo de area de nodo (ha)
# Fundamento: Rabinowitz & Zeller (2010); Noss et al. (2003)
UMBRAL_HA = 500.0

RES = 90.0                      # resolucion (m)
AP  = RES ** 2                  # area pixel (m2)
AH  = AP / 10_000               # area pixel (ha)
S8  = ndimage.generate_binary_structure(2, 2)   # 8-conectividad


def extraer_nodos(anio: int):
    """
    Extrae nodos nucleo >= UMBRAL_HA para un anio y genera raster ASC.

    Args:
        anio: Anio de analisis.

    Returns:
        Dict con id_nodo -> {area_ha, centroide_xy, bbox}.
    """
    ruta_in   = os.path.join(IN_DIR, f"mb_riogrande_{anio}_90m.tif")
    ruta_asc  = os.path.join(OUT_DIR, f"nodos_{anio}_90m.asc")
    ruta_csv  = os.path.join(CSV_DIR, f"nodos_atributos_{anio}.csv")

    ds   = gdal.Open(ruta_in)
    if ds is None:
        raise FileNotFoundError(f"No se encontro: {ruta_in}")
    arr  = ds.GetRasterBand(1).ReadAsArray().astype(np.int16)
    nd   = int(ds.GetRasterBand(1).GetNoDataValue() or -9999)
    gt   = ds.GetGeoTransform()
    proj = ds.GetProjection()
    rows, cols = arr.shape
    ds = None

    # Mascara de bosque
    bin_a = np.zeros_like(arr, dtype=np.uint8)
    for c in CODIGOS_BOSQUE:
        bin_a[arr == c] = 1
    bin_a[arr == nd] = 0

    # Etiquetar parches
    lab, n_total = ndimage.label(bin_a, structure=S8)
    sz_ha = np.array(ndimage.sum(bin_a, lab, range(1, n_total + 1))) * AH

    # Filtrar por umbral
    ids_validos = [i + 1 for i, s in enumerate(sz_ha) if s >= UMBRAL_HA]
    print(f"  {anio}: {n_total} parches totales | {len(ids_validos)} nodos >= {UMBRAL_HA} ha")

    # Raster de nodos (ID incremental)
    nodo_arr   = np.zeros_like(arr, dtype=np.int16)
    atributos  = {}
    for nuevo_id, pid in enumerate(sorted(ids_validos,
                                          key=lambda x: sz_ha[x-1],
                                          reverse=True), start=1):
        mascara = lab == pid
        nodo_arr[mascara] = nuevo_id
        cy, cx  = ndimage.center_of_mass(mascara)
        # Convertir pixeles a coordenadas geograficas
        x_geo   = gt[0] + cx * gt[1]
        y_geo   = gt[3] + cy * gt[5]
        area_ha = round(float(sz_ha[pid - 1]), 2)
        atributos[nuevo_id] = {
            "nodo_id"  : nuevo_id,
            "pid_orig" : pid,
            "area_ha"  : area_ha,
            "centX"    : round(x_geo, 1),
            "centY"    : round(y_geo, 1),
            "anio"     : anio,
        }

    # Exportar CSV de atributos
    if atributos:
        with open(ruta_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(next(iter(atributos.values())).keys()))
            writer.writeheader()
            writer.writerows(atributos.values())

    # Exportar raster ASC para Circuitscape
    # ASC = formato ASCII Grid requerido por Circuitscape 4.0
    xmin = gt[0]
    ymax = gt[3]
    dx   = gt[1]
    with open(ruta_asc, "w") as f:
        f.write(f"ncols         {cols}\n")
        f.write(f"nrows         {rows}\n")
        f.write(f"xllcorner     {xmin:.6f}\n")
        f.write(f"yllcorner     {ymax + rows * gt[5]:.6f}\n")
        f.write(f"cellsize      {dx:.2f}\n")
        f.write(f"NODATA_value  -9999\n")
        # Escribir por bloques para eficiencia
        for row in nodo_arr:
            f.write(" ".join(str(v) for v in row) + "\n")

    print(f"         -> ASC: {os.path.basename(ruta_asc)}")
    print(f"         -> CSV: {os.path.basename(ruta_csv)}")
    return atributos


# ── Ejecutar ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    ANIOS_CS = [1985, 2005, 2024]

    print("=" * 65)
    print("EXTRACCION DE NODOS NUCLEO — Subcuenca rio Grande-Pirai")
    print(f"Umbral: >= {UMBRAL_HA} ha | Conectividad: 8 celdas")
    print("=" * 65)

    for anio in ANIOS_CS:
        extraer_nodos(anio)

    print("\nOK: nodos extraidos. Siguiente paso:")
    print("    Ejecutar Circuitscape con los .ini de /circuitscape/")
    print("    o ejecutar directamente 04_conectividad.py")
