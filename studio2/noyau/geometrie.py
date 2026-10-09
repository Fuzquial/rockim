"""Géométrie des discontinuités et estimations de maillage, sans dépendance."""
import math


def segments_plans(pendage_deg, espacement, fraction_rompue, W, H):
    """Segments (x1, y1, x2, y2) des plans DÉJÀ ROMPUS d'une famille.

    Repris tel quel de etude_triax_hetero/gen_decks.py (plane_segments), qui suit
    la convention de FdemSolver::applyWeakPlanes : normale n = (-sin b, cos b),
    plans n . X = k * espacement ANCRÉS SUR L'ORIGINE du maillage. Ancrés au
    centre, les segments tombaient à mi-chemin des plans affaiblis (mesure du
    2026-09-15). Les cordes de moins de 10 % de la diagonale sont écartées ; la
    fraction rompue est prise régulièrement répartie dans la famille.
    """
    b = math.radians(pendage_deg)
    nx, ny = -math.sin(b), math.cos(b)
    kmax = int(math.ceil((abs(nx) * W + abs(ny) * H) / espacement)) + 2
    segs = []
    for k in range(-kmax, kmax + 1):
        d = k * espacement
        pts = []
        for x in (0.0, W):
            if abs(ny) > 1e-12:
                y = (d - nx * x) / ny
                if -1e-12 <= y <= H + 1e-12:
                    pts.append((x, min(max(y, 0.0), H)))
        for y in (0.0, H):
            if abs(nx) > 1e-12:
                x = (d - ny * y) / nx
                if -1e-12 <= x <= W + 1e-12:
                    pts.append((min(max(x, 0.0), W), y))
        if len(pts) < 2:
            continue
        (x1, y1), (x2, y2) = pts[0], pts[-1]
        if math.hypot(x2 - x1, y2 - y1) < 0.10 * math.hypot(W, H):
            continue
        segs.append((x1, y1, x2, y2))
    if not segs or fraction_rompue <= 0.0:
        return []
    n = max(1, int(round(fraction_rompue * len(segs))))
    if n >= len(segs):
        return segs
    pas = len(segs) / float(n)
    return [segs[min(len(segs) - 1, int(round((i + 0.5) * pas)))] for i in range(n)]


def elements_par_grain(taille_grain, taille_element):
    """Estimation du nombre de triangles par grain de Voronoï.

    Aire moyenne d'un grain ~ (pi / 4) d², aire d'un triangle de Delaunay de
    côté h ~ (sqrt 3 / 4) h², d'où N ~ (pi / sqrt 3) (d / h)². Contrôle sur la
    campagne (d = 3 mm, h = 0,7 mm, log de F2_gbm_P020) : 367 grains et 12 995
    éléments sur 36 x 72 mm, soit 35,4 mesurés pour 33,3 estimés (-6 %).
    On retient le coefficient MESURÉ, 1,93 au lieu de pi / sqrt 3 = 1,81 :
    l'estimation sert à appliquer la règle maison de 35-90 éléments par grain
    (2026-09-07), et la formule brute ferait échouer la campagne qui la respecte.
    """
    return K_ELEMENTS_PAR_GRAIN * (taille_grain / taille_element) ** 2


# 12 995 éléments / 367 grains / (3 / 0,7)² ; log de F2_gbm_P020 (2026-09-15).
K_ELEMENTS_PAR_GRAIN = 12995.0 / 367.0 / (0.003 / 0.0007) ** 2


def aire(essai):
    """Aire de l'éprouvette : rectangle W x H, ou disque à méplats du brésilien. Chaque méplat
    retire un segment circulaire d'angle theta = discFlattenDeg (angle TOTAL) :
    A = R^2 (pi - theta + sin theta)."""
    e, c = essai.eprouvette, essai.chargement
    if c.type_essai == "bresilien":
        R, th = e.W / 2, math.radians(c.aplatissement_deg)
        return R * R * (math.pi - th + math.sin(th))
    return e.W * e.H


def nombre_elements(W, H, taille_element):
    """Triangles de Delaunay de côté h sur W x H. Contrôle F2 : 12 995 mesurés,
    W H / (0,41 h²) = 12 990 avec le coefficient 0,407 mesuré (aire moyenne)."""
    return W * H / (0.407 * taille_element ** 2)
