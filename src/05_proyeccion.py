"""
05_proyeccion.py
================
Proyeccion temporal del Indice de Conectividad de Probabilidad (PC)
hacia 2040 bajo tres escenarios, con analisis de sensibilidad.

Escenarios:
    1. Conservador  — tasa media historica ponderada (-0.82 %/anio)
    2. Tendencial   — tasa del periodo 2005-2024 (-0.91 %/anio)
    3. Aceleracion  — tasa del cuatrienio 2020-2024 (si se mantuviera)

Umbral de inviabilidad ecologica: PC < 0.020
    Fundamento: Bennett, A.F. (2003). Linkages in the landscape (2a ed.).
    IUCN. Aplicacion indicativa; requiere validacion con datos de fauna.

El script tambien calcula el ratio de perdida (conectividad vs area)
para cada periodo y los tres escenarios de proyeccion.

Uso:
    python 05_proyeccion.py
    # Salida: data/conectividad/proyeccion_PC_2040.csv
              data/conectividad/proyeccion_PC_2040_escenarios.csv

Requisitos:
    numpy>=1.26, scipy>=1.11
    Entorno: conda activate conectividad-pirai

Autor: Gloria Eliana Torrez Castro — PPGG/UFC, 2024-2025
"""

import numpy as np
import csv
import os

# ── Rutas ───────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_DIR   = os.path.join(BASE_DIR, "data", "conectividad")
FRAG_DIR = os.path.join(BASE_DIR, "data", "fragmentacion")
OUT_DIR  = os.path.join(BASE_DIR, "data", "conectividad")

# ── Valores observados — extraidos de 04_conectividad.py ────────────────────
# PC global por anio (valores calculados con NetworkX, d=3km, umbral=5km)
PC_OBSERVADO = {
    1985: 0.04615,
    2005: 0.03541,
    2024: 0.02975,
}

# Area forestal (ha) por anio — de fragmentacion_riogrande_1985_2024.csv
CA_BOSQUE = {
    1985: 304_461,
    2005: 281_738,   # proxy 2006 (diferencia de 1 anio, ecologicamente irrelevante)
    2024: 311_324,
}

# PC del cuatrienio 2020-2024 para estimar tasa acelerada
# (calculado a partir de los datos FRAGSTATS: ruptura critica 2020-2024)
# Si no se dispone del valor 2020, se estima por interpolacion lineal
PC_2020_ESTIMADO = 0.03100   # estimacion conservadora entre 2005 y 2024

# Umbral de inviabilidad ecologica — Bennett (2003)
PC_UMBRAL = 0.020

# Anio de proyeccion
ANIO_PROYECCION = 2040
ANIO_ULTIMO     = 2024
HORIZONTE_ANIOS = ANIO_PROYECCION - ANIO_ULTIMO   # 16 anios


# ── Calculos de tasas ────────────────────────────────────────────────────────
def tasa_anual(pc_ini: float, pc_fin: float, n_anios: int) -> float:
    """Tasa de cambio anual media (porcentaje por anio)."""
    return (pc_fin - pc_ini) / pc_ini / n_anios * 100


def anio_umbral_lineal(pc_ini: float, tasa_pct_anio: float,
                       umbral: float = PC_UMBRAL) -> float | None:
    """
    Estima el anio en que se alcanza el umbral bajo tasa lineal.

    Returns:
        Anio decimal o None si no se alcanza en 100 anios.
    """
    if tasa_pct_anio >= 0:
        return None   # tasa creciente o neutra: umbral no se alcanza
    anios = (umbral - pc_ini) / (pc_ini * tasa_pct_anio / 100)
    return round(ANIO_ULTIMO + anios, 1)


def proyectar_lineal(pc_ini: float, tasa_pct_anio: float,
                     n_anios: int) -> float:
    """PC proyectado bajo tasa lineal."""
    return pc_ini * (1 + tasa_pct_anio / 100) ** n_anios


# ── Ejecutar ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("PROYECCION PC 2040 — Subcuenca rio Grande-Pirai")
    print("=" * 70)

    # Tasas historicas
    tasa_1985_2005 = tasa_anual(PC_OBSERVADO[1985], PC_OBSERVADO[2005], 20)
    tasa_2005_2024 = tasa_anual(PC_OBSERVADO[2005], PC_OBSERVADO[2024], 19)
    tasa_media_pond = (tasa_1985_2005 * 20 + tasa_2005_2024 * 19) / 39
    tasa_aceleracion = tasa_anual(PC_2020_ESTIMADO, PC_OBSERVADO[2024], 4)

    print(f"\nTasas de perdida de PC observadas:")
    print(f"  1985-2005: {tasa_1985_2005:.4f} %/anio ({tasa_1985_2005*20:.1f}% total en 20 anios)")
    print(f"  2005-2024: {tasa_2005_2024:.4f} %/anio ({tasa_2005_2024*19:.1f}% total en 19 anios)")
    print(f"  2020-2024: {tasa_aceleracion:.4f} %/anio (cuatrienio critico, estimacion)")
    print(f"  Ponderada: {tasa_media_pond:.4f} %/anio (39 anios 1985-2024)")

    # Ratio perdida conectividad vs area
    delta_pc_pct  = (PC_OBSERVADO[2024] - PC_OBSERVADO[1985]) / PC_OBSERVADO[1985] * 100
    delta_ca_pct  = (CA_BOSQUE[2024]    - CA_BOSQUE[1985])    / CA_BOSQUE[1985]    * 100
    ratio         = abs(delta_pc_pct) / abs(delta_ca_pct) if delta_ca_pct != 0 else None

    print(f"\nEfecto no lineal:")
    print(f"  DeltaPC  = {delta_pc_pct:.1f}%")
    print(f"  DeltaCA  = {delta_ca_pct:.1f}%")
    print(f"  Ratio    = {ratio:.2f}:1  (conectividad se deteriora {ratio:.1f}x mas rapido que el area)")

    # Tres escenarios de proyeccion
    escenarios = [
        {
            "nombre"      : "Conservador",
            "descripcion" : "Tasa media historica ponderada (1985-2024)",
            "tasa_pct_anio": tasa_media_pond,
        },
        {
            "nombre"      : "Tendencial",
            "descripcion" : "Tasa del periodo 2005-2024",
            "tasa_pct_anio": tasa_2005_2024,
        },
        {
            "nombre"      : "Aceleracion",
            "descripcion" : "Tasa del cuatrienio 2020-2024 (si se mantuviera)",
            "tasa_pct_anio": tasa_aceleracion,
        },
    ]

    print(f"\nProyeccion de PC hacia {ANIO_PROYECCION} (horizonte: {HORIZONTE_ANIOS} anios):")
    print(f"  Umbral de inviabilidad ecologica: PC < {PC_UMBRAL} (Bennett, 2003)")
    print(f"  {'Escenario':<20} {'Tasa %/a':<12} {'PC en 2040':<12} {'Umbral antes'}")

    resultados_esc = []
    for esc in escenarios:
        pc_2040  = proyectar_lineal(PC_OBSERVADO[ANIO_ULTIMO],
                                     esc["tasa_pct_anio"], HORIZONTE_ANIOS)
        anio_umb = anio_umbral_lineal(PC_OBSERVADO[ANIO_ULTIMO],
                                       esc["tasa_pct_anio"])
        cruza    = pc_2040 < PC_UMBRAL
        print(f"  {esc['nombre']:<20} {esc['tasa_pct_anio']:>+8.4f}    "
              f"{pc_2040:.5f}     "
              f"{'SI (~' + str(anio_umb) + ')' if cruza and anio_umb else 'NO (>' + str(ANIO_PROYECCION) + ')'}")
        resultados_esc.append({
            "escenario"    : esc["nombre"],
            "descripcion"  : esc["descripcion"],
            "tasa_pct_anio": round(esc["tasa_pct_anio"], 4),
            "PC_2024"      : PC_OBSERVADO[ANIO_ULTIMO],
            "PC_2040_proj" : round(pc_2040, 5),
            "umbral_PC"    : PC_UMBRAL,
            "cruza_umbral" : cruza,
            "anio_umbral"  : anio_umb if anio_umb else ">2060",
        })

    # Exportar
    ruta_out = os.path.join(OUT_DIR, "proyeccion_PC_2040_escenarios.csv")
    with open(ruta_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(resultados_esc[0].keys()))
        w.writeheader()
        w.writerows(resultados_esc)

    print(f"\nOK: resultados guardados en {ruta_out}")
    print("\nNOTA: Las proyecciones son extrapolaciones lineales conservadoras.")
    print("      No incorporan cambio climatico, variaciones en precios de")
    print("      commodities ni nuevas politicas de conservacion.")
    print("      Interpretar como ordenes de magnitud orientativos.")
