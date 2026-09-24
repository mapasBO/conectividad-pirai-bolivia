"""
01_fragmentacion.py
===================
Metricas de fragmentacion estructural del bosque — Subcuenca rio Grande-Pirai
Equivalente a FRAGSTATS (McGarigal & Marks, 1995) con SciPy/NumPy.

Metricas calculadas por clase y anio:
    NP    — Numero de parches
    CA    — Area de clase (ha)
    LPI   — Largest Patch Index (%)
    MPS   — Mean Patch Size (ha)
    PSSD  — Patch Size Standard Deviation (ha)
    FRAC  — Mean Fractal Dimension
    PARA  — Perimeter-Area Ratio medio
    ENN   — Mean Euclidean Nearest-Neighbor Distance (km)
    SHDI  — Shannon Diversity Index del paisaje completo

Referencia: McGarigal, K. & Marks, B.J. (1995). FRAGSTATS. USDA-FS GTR PNW-351.

Uso:
    python 01_fragmentacion.py
    # o desde Consola Python de QGIS: exec(open('01_fragmentacion.py').read())

Requisitos:
    numpy>=1.26, scipy>=1.11, gdal>=3.12
    Entorno: conda activate conectividad-pirai  (ver environment.yml)

Autor: Gloria Eliana Torrez Castro — PPGG/UFC, 2024–2025
"""

import numpy as np
import os
import csv
import time
from scipy import ndimage
from scipy.spatial import cKDTree
from osgeo import gdal

# ── Rutas ───────────────────────────────────────────────────────────────────
# Adaptar BASE_DIR a la instalacion local
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_DIR   = os.path.join(BASE_DIR, "data", "rasters_30m")   # rasters MapBiomas recortados
OUT_DIR  = os.path.join(BASE_DIR, "data", "fragmentacion")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parametros ───────────────────────────────────────────────────────────────
ANIOS           = [1985, 2001, 2005, 2010, 2015, 2020, 2024]
RES             = 30.0          # resolucion espacial (m)
AP              = RES ** 2      # area pixel (m2)
AH              = AP / 10_000   # area pixel (ha)
MAX_ENN_SAMPLE  = 300           # subsample para ENN (eficiencia computacional)
MAX_FRAC_SAMPLE = 150           # top parches por area para FRAC
S8              = ndimage.generate_binary_structure(2, 2)   # 8-conectividad

# Clases MapBiomas Bolivia Coleccion 2024 — codigos Nivel 2-3
# Fuente: https://bolivia.mapbiomas.org/codigos-de-la-leyenda/
CLASES = {
    "Bosque"         : [3, 6, 11, 68, 72],   # Formacion bocosa + humedal forestal
    "Agropecuario"   : [15, 18, 21],
    "Herbazal"       : [12, 13],
    "Sin_vegetacion" : [23, 24, 25, 30],
    "Agua"           : [31, 33],
}

# Umbral ecologico de referencia MPS (ha) — Bennett (2003)
UMBRAL_MPS_HA = 50.0


# ── Funcion principal ────────────────────────────────────────────────────────
def analizar_anio(anio: int) -> list:
    """
    Calcula metricas FRAGSTATS para todas las clases en un anio.

    Args:
        anio: Anio de analisis (debe existir mb_riogrande_{anio}.tif en IN_DIR)

    Returns:
        Lista de dicts con metricas por clase.
    """
    t0   = time.time()
    ruta = os.path.join(IN_DIR, f"mb_riogrande_{anio}.tif")
    ds   = gdal.Open(ruta)
    if ds is None:
        raise FileNotFoundError(f"No se encontro: {ruta}")
    arr  = ds.GetRasterBand(1).ReadAsArray().astype(np.int16)
    nd   = int(ds.GetRasterBand(1).GetNoDataValue() or -9999)
    ds   = None

    msk = arr != nd
    tp  = int(msk.sum())           # total pixeles validos
    am  = np.where(msk, arr, 0)
    th  = tp * AH                  # area total valida (ha)

    # SHDI — indice de diversidad Shannon del paisaje completo
    props = []
    for cods in CLASES.values():
        px = int(sum((am == c).sum() for c in cods))
        if px > 0:
            props.append(px / tp)
    SHDI = round(-sum(p * np.log(p) for p in props if p > 0), 4)

    filas = []
    for cls, cods in CLASES.items():
        bin_a = np.zeros_like(am, dtype=np.uint8)
        for c in cods:
            bin_a[am == c] = 1
        px_cls = int(bin_a.sum())
        if px_cls == 0:
            continue

        CA     = round(px_cls * AH, 1)
        CA_pct = round(CA / th * 100, 3)

        lab, NP = ndimage.label(bin_a, structure=S8)
        sz_ha   = np.array(ndimage.sum(bin_a, lab, range(1, NP + 1))) * AH
        LPI     = round(float(np.max(sz_ha)) / th * 100, 4)
        MPS     = round(float(np.mean(sz_ha)), 3)
        PSSD    = round(float(np.std(sz_ha)), 1)

        # ENN — distancia Euclidiana media al vecino mas proximo (km)
        ids_s = (range(1, NP + 1) if NP <= MAX_ENN_SAMPLE
                 else list(np.random.choice(range(1, NP + 1),
                                            MAX_ENN_SAMPLE, replace=False)))
        cents = np.array([ndimage.center_of_mass(bin_a, lab, i)
                          for i in ids_s]) * RES
        if len(cents) > 1:
            d2, _   = cKDTree(cents).query(cents, k=2)
            ENN_mn  = round(float(np.mean(d2[:, 1])) / 1000, 4)
            ENN_cv  = round((float(np.std(d2[:, 1])) / float(np.mean(d2[:, 1]))
                             if np.mean(d2[:, 1]) > 0 else 0), 4)
        else:
            ENN_mn = ENN_cv = 0.0

        # FRAC — dimension fractal media (top parches por area)
        top = list(np.argsort(sz_ha)[-MAX_FRAC_SAMPLE:] + 1)
        fracs, paras = [], []
        for pid in top:
            p2  = lab == pid
            a2  = float(ndimage.sum(p2)) * AP
            pm  = float((ndimage.binary_dilation(p2) & ~p2).sum()) * RES
            if a2 > AP and pm > 0:
                f = 2 * np.log(pm / 4) / np.log(a2)
                fracs.append(min(max(f, 1.0), 2.0))
                paras.append(pm / a2)
        FRAC = round(float(np.mean(fracs)) if fracs else 1.0, 4)
        PARA = round(float(np.mean(paras)) if paras else 0.0, 6)

        filas.append({
            "anio"      : anio,
            "clase"     : cls,
            "NP"        : NP,
            "CA_ha"     : CA,
            "CA_pct"    : CA_pct,
            "LPI_pct"   : LPI,
            "MPS_ha"    : MPS,
            "PSSD"      : PSSD,
            "FRAC_mn"   : FRAC,
            "PARA_mn"   : PARA,
            "ENN_mn_km" : ENN_mn,
            "ENN_cv"    : ENN_cv,
            "SHDI"      : SHDI,
            "total_ha"  : round(th, 0),
            "umbral_MPS": "BAJO UMBRAL" if (cls == "Bosque"
                                             and MPS < UMBRAL_MPS_HA) else "OK",
        })

    t1 = time.time()
    b  = next((f for f in filas if f["clase"] == "Bosque"), {})
    print(f"  {anio} ({t1-t0:.0f}s) | Bosque: NP={b.get('NP',0):>5,} "
          f"CA={b.get('CA_ha',0):>9,.0f} ha  MPS={b.get('MPS_ha',0):.1f} ha "
          f"ENN={b.get('ENN_mn_km',0):.3f} km  FRAC={b.get('FRAC_mn',0):.4f}  "
          f"SHDI={SHDI:.4f}  {b.get('umbral_MPS','')}")
    return filas


# ── Ejecutar ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 75)
    print("ANALISIS DE FRAGMENTACION — Subcuenca rio Grande-Pirai")
    print("MapBiomas Bolivia Coleccion 2024 | 7 cortes temporales")
    print("=" * 75)

    todas = []
    for anio in ANIOS:
        print(f"  Procesando {anio}...")
        todas.extend(analizar_anio(anio))

    ruta_csv = os.path.join(OUT_DIR, "fragmentacion_riogrande_1985_2024.csv")
    campos   = list(todas[0].keys())
    with open(ruta_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(todas)

    print(f"\n{'='*75}")
    print(f"OK: {len(todas)} registros exportados")
    print(f"OK: CSV guardado en {ruta_csv}")
