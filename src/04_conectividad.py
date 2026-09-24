"""
04_conectividad.py
==================
Calculo de indices de conectividad funcional con NetworkX.
Equivalente a Conefor Sensinode 2.2 (Saura & Torne, 2009).

Indices calculados:
    IIC  — Integral Index of Connectivity
           (Pascual-Hortal & Saura, 2006)
    PC   — Probability of Connectivity
           (Saura & Pascual-Hortal, 2007)
    dPC  — Variacion porcentual de PC al eliminar cada nodo
           (indicador de importancia nodal para priorizacion)

Parametros del grafo:
    - Umbral de conexion directa: 5 km (distancia Euclidiana entre centroides)
    - Funcion de dispersion: p(d) = exp(-1 * d / d_media)
      con d_media = 3 km (distancia de dispersion efectiva para tapir y jaguar)
    - Area del paisaje: A_L = 1_352_700 ha (subcuenca rio Grande-Pirai)

Referencia metodologica:
    Saura, S. & Torne, J. (2009). Conefor Sensinode 2.2.
    Env. Modelling & Software 24(1): 135-139.
    https://doi.org/10.1016/j.envsoft.2008.05.005

    Hagberg, A., Swart, P. & Chult, D. (2008). Exploring network structure,
    dynamics, and function using NetworkX (LA-UR-08-05495). LANL.

Uso:
    python 04_conectividad.py
    # Salida: data/conectividad/indices_PC_IIC_dPC.csv

Requisitos:
    numpy>=1.26, networkx>=3.2, scipy>=1.11
    Entorno: conda activate conectividad-pirai

Autor: Gloria Eliana Torrez Castro — PPGG/UFC, 2024-2025
"""

import numpy as np
import networkx as nx
import csv
import os
import math
from itertools import combinations

# ── Rutas ───────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_DIR  = os.path.join(BASE_DIR, "data", "fragmentacion")
OUT_DIR  = os.path.join(BASE_DIR, "data", "conectividad")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parametros del modelo de conectividad ────────────────────────────────────
D_MEDIA_KM      = 3.0          # distancia de dispersion efectiva (km)
                                # Fundamento: home range tapir 200-500 ha en
                                # paisajes fragmentados (Noss et al., 2003);
                                # distancia minima 2 km, maxima 5 km.
UMBRAL_CONN_KM  = 5.0          # umbral de conexion directa (km)
                                # Parches separados por > 5 km no se conectan
                                # directamente en el grafo.
AL_HA           = 1_352_700.0  # area del paisaje (ha) — subcuenca completa

# Anos de analisis de grafos (subset de los 7 cortes FRAGSTATS)
# Fundamento de seleccion: representatividad de periodos historico-politicos
# contrastantes (neoliberal 1985-2005, plurinacional 2005-2024) y reduccion
# de carga computacional del calculo iterativo de dPC.
ANIOS_GRAFOS = [1985, 2005, 2024]


# ── Funciones de conectividad ────────────────────────────────────────────────
def prob_dispersion(dist_km: float, d_media: float = D_MEDIA_KM) -> float:
    """
    Probabilidad de dispersion con decaimiento exponencial negativo.
    p(d) = exp(-1 * d / d_media)

    Args:
        dist_km: Distancia entre centroides de parches (km).
        d_media: Distancia media de dispersion efectiva (km).

    Returns:
        Probabilidad [0, 1].
    """
    return math.exp(-dist_km / d_media)


def construir_grafo(atributos: dict) -> nx.Graph:
    """
    Construye el grafo de conectividad a partir de los atributos de nodos.

    Args:
        atributos: Dict {nodo_id: {area_ha, centX, centY, ...}}

    Returns:
        Grafo no dirigido ponderado con pesos = probabilidad de dispersion.
    """
    G = nx.Graph()
    for nid, attrs in atributos.items():
        G.add_node(nid, area_ha=attrs["area_ha"],
                   x=attrs["centX"], y=attrs["centY"])

    nodos = list(atributos.keys())
    for n1, n2 in combinations(nodos, 2):
        x1, y1 = atributos[n1]["centX"], atributos[n1]["centY"]
        x2, y2 = atributos[n2]["centX"], atributos[n2]["centY"]
        dist_m  = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        dist_km = dist_m / 1000

        if dist_km <= UMBRAL_CONN_KM:
            p = prob_dispersion(dist_km)
            G.add_edge(n1, n2, dist_km=round(dist_km, 3), prob=round(p, 6))

    return G


def calcular_IIC(G: nx.Graph, AL: float) -> float:
    """
    Calcula el Integral Index of Connectivity (IIC).
    IIC = sum_i sum_j (ai * aj) / (1 + nlinks_ij) / AL^2

    Args:
        G:  Grafo de conectividad.
        AL: Area del paisaje (ha).

    Returns:
        Valor de IIC [0, 1].
    """
    nodos   = list(G.nodes)
    total   = 0.0
    paths   = dict(nx.all_pairs_shortest_path_length(G))
    for i in nodos:
        ai = G.nodes[i]["area_ha"]
        for j in nodos:
            aj = G.nodes[j]["area_ha"]
            if j in paths.get(i, {}):
                nlinks = paths[i][j]
            else:
                nlinks = 1e9   # sin conexion -> contribucion ~0
            total += (ai * aj) / (1 + nlinks)
    return total / (AL ** 2)


def calcular_PC(G: nx.Graph, AL: float) -> float:
    """
    Calcula el Probability of Connectivity (PC).
    PC = sum_i sum_j (ai * aj * pij*) / AL^2
    donde pij* = maxima probabilidad de dispersion entre i y j
    sobre todos los caminos posibles.

    Args:
        G:  Grafo de conectividad (pesos = probabilidad directa).
        AL: Area del paisaje (ha).

    Returns:
        Valor de PC [0, 1].
    """
    nodos = list(G.nodes)
    # Construir grafo con log de probabilidades para usar Dijkstra
    G_log = nx.Graph()
    G_log.add_nodes_from(G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        # -log(p) convierte el problema a minimizacion de camino minimo
        w = -math.log(d["prob"]) if d["prob"] > 0 else 1e9
        G_log.add_edge(u, v, weight=w)

    total = 0.0
    for i in nodos:
        ai = G.nodes[i]["area_ha"]
        lengths = nx.single_source_dijkstra_path_length(G_log, i, weight="weight")
        for j in nodos:
            aj = G.nodes[j]["area_ha"]
            if j in lengths:
                pij = math.exp(-lengths[j])
            else:
                pij = 0.0
            total += ai * aj * pij
    return total / (AL ** 2)


def calcular_dPC(G: nx.Graph, AL: float, pc_global: float) -> dict:
    """
    Calcula la importancia nodal dPC para cada nodo.
    dPC_i = (PC_global - PC_sin_i) / PC_global * 100

    Args:
        G:         Grafo completo.
        AL:        Area del paisaje (ha).
        pc_global: PC calculado con todos los nodos.

    Returns:
        Dict {nodo_id: dPC_porcentaje}
    """
    nodos = list(G.nodes)
    resultado = {}
    for i, nodo in enumerate(nodos):
        G_sin = G.copy()
        G_sin.remove_node(nodo)
        pc_sin = calcular_PC(G_sin, AL) if len(G_sin.nodes) > 0 else 0.0
        dpc    = (pc_global - pc_sin) / pc_global * 100 if pc_global > 0 else 0.0
        resultado[nodo] = round(dpc, 4)
        print(f"    dPC N{nodo:02d}: {dpc:.4f}% [{i+1}/{len(nodos)}]")
    return resultado


def leer_atributos_csv(anio: int) -> dict:
    """Lee el CSV de atributos de nodos generado por 03_nodos.py."""
    ruta = os.path.join(CSV_DIR, f"nodos_atributos_{anio}.csv")
    atributos = {}
    with open(ruta, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            nid = int(row["nodo_id"])
            atributos[nid] = {
                "area_ha": float(row["area_ha"]),
                "centX"  : float(row["centX"]),
                "centY"  : float(row["centY"]),
            }
    return atributos


# ── Ejecutar ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("INDICES DE CONECTIVIDAD FUNCIONAL — IIC, PC, dPC")
    print("Subcuenca rio Grande-Pirai | NetworkX 3.2")
    print(f"Parametros: d_media={D_MEDIA_KM} km | umbral={UMBRAL_CONN_KM} km | AL={AL_HA:,.0f} ha")
    print("=" * 70)

    resultados = []

    for anio in ANIOS_GRAFOS:
        print(f"\n--- {anio} ---")
        atributos = leer_atributos_csv(anio)
        n_nodos   = len(atributos)
        print(f"  Nodos: {n_nodos} | Iteraciones dPC: {n_nodos}")

        G   = construir_grafo(atributos)
        iic = calcular_IIC(G, AL_HA)
        pc  = calcular_PC(G, AL_HA)
        print(f"  IIC = {iic:.5f} | PC = {pc:.5f}")

        print("  Calculando dPC (puede tardar 2-5 minutos por anio)...")
        dpc_dict = calcular_dPC(G, AL_HA, pc)

        for nod, attrs in atributos.items():
            resultados.append({
                "anio"   : anio,
                "nodo_id": nod,
                "area_ha": attrs["area_ha"],
                "centX"  : attrs["centX"],
                "centY"  : attrs["centY"],
                "n_edges": G.degree(nod),
                "IIC"    : round(iic, 6),
                "PC"     : round(pc, 6),
                "dPC_pct": dpc_dict.get(nod, 0.0),
                "AL_ha"  : AL_HA,
            })

    # Exportar CSV
    ruta_out = os.path.join(OUT_DIR, "indices_PC_IIC_dPC.csv")
    campos   = list(resultados[0].keys())
    with open(ruta_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(resultados)

    print(f"\n{'='*70}")
    print(f"OK: {len(resultados)} registros exportados")
    print(f"OK: CSV guardado en {ruta_out}")
    print("\nSiguiente paso: ejecutar 05_proyeccion.py")
