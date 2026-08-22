# ---------------------------------------------------------------------------
# modules.py — les MODULES métier du studio (demande F. Uzquiano 2026-08-22).
#
# Au démarrage on choisit un métier (tunnel EDZ, hydro-frac, impact, essais
# de laboratoire) et l'interface ne montre plus QUE ce qui lui appartient :
# ses gabarits, ses groupes de clés (le panneau Matériau en tête), et ses
# actions propres — le balayage lambda pour le tunnel. Le mode « Expert »
# garde l'interface complète (comportement historique).
#
# Tout est DÉCLARATIF : un module = gabarits + groupes autorisés + balayages.
# Un groupe hors liste n'est affiché que s'il porte des clés posées (on ne
# cache jamais une clé explicitement présente dans le cfg).
# ---------------------------------------------------------------------------
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Sweep:
    """Un balayage paramétrique à un facteur, façon Wang et al. 2024."""
    name: str                    # étiquette du menu
    doc: str                     # explication affichée dans le dialogue
    values: str                  # valeurs proposées par défaut ("0.5 0.75 …")
    # make(cas, v) -> dict des clés à poser pour la valeur v
    make: object = None


def _lambda_keys(v: float) -> dict[str, str]:
    """Convention validée sur leur fig. 16 (make_configs.py du 2026-08-17) :
    sigma_h = 5 MPa FIXE, sigma_v = 5/lambda — le nombre de fissures décroît
    avec lambda comme chez eux."""
    return {"insituSh": "5e6", "insituSv": "%.10g" % (5e6 / v)}


def _sigma0_keys(v: float) -> dict[str, str]:
    return {"insituSh": "%.10g" % (v * 1e6), "insituSv": "%.10g" % (v * 1e6)}


@dataclass(frozen=True)
class Module:
    key: str
    label: str
    doc: str
    templates: tuple = ()
    groups: tuple = ()           # () = tous (mode expert)
    sweeps: tuple = ()


MODULES = (
    Module(
        "tunnel", "Tunnel EDZ",
        "Excavation, contraintes in situ, zone endommagée "
        "(Wang et al. 2024 — Hutou Beishan).",
        templates=(
            ("Tunnel Wang 2024 (référence σ₀ = 5 MPa, λ = 1)",
             "tunnel_edz/configs/tunnel_ref_s5_lam1.cfg"),
            ("Tunnel EDZ pressurisé (rapide)", "configs/tunnel_bore_fast.cfg"),
            ("Tunnel EDZ pressurisé (production)", "configs/tunnel_bore.cfg"),
            ("Tunnel EDZ Weibull", "configs/tunnel_bore_weib.cfg"),
        ),
        groups=("Général", "Maillage", "Matériau", "Joints", "Contact",
                "Conditions aux limites", "Sorties"),
        sweeps=(
            Sweep("Balayage λ (coefficient de pression latérale)…",
                  "Leur fig. 15-16 : σ_h = 5 MPa fixe, σ_v = 5/λ. "
                  "Un deck et un run par valeur, mis en file l'un "
                  "après l'autre.",
                  "0.5 0.75 1.25 1.5", _lambda_keys),
            Sweep("Balayage σ₀ (in situ hydrostatique)…",
                  "Leur fig. 12-14 : σ_h = σ_v = σ₀ [MPa].",
                  "3 4 6 7", _sigma0_keys),
        ),
    ),
    Module(
        "hydro", "Fracturation hydraulique",
        "Pompe, cavité, pression de rupture (AbuAisha 2017, spec 004).",
        templates=(
            ("ISO grossier (hydro = on, ~6 min)",
             "bench_abuaisha/configs/hf_iso_hydro_c.cfg"),
            ("ISO production", "bench_abuaisha/configs/hf_iso_hydro_s.cfg"),
            ("ANISO production", "bench_abuaisha/configs/hf_aniso_hydro_s.cfg"),
            ("Essai 3 : protocole article (hydroStart)",
             "bench_abuaisha/configs/e3_iso12.cfg"),
        ),
        groups=("Général", "Maillage", "Matériau", "Joints", "Contact",
                "Conditions aux limites", "Hydro", "Sorties"),
    ),
    Module(
        "impact", "Impact / percussion",
        "Frappe à insert, multi-corps, pulvérisation "
        "(Yang et al. 2025-26, spec 005).",
        templates=(
            ("Impact 3D smoke (bench1 réduit)", "configs/smoke_impact.cfg"),
            ("Impact banc St Anne s1,5 (spec 005)",
             "bench_impact/configs/impact_stanne_s15.cfg"),
            ("Impact St Anne FIDÈLE (tout comme eux sauf adaptatif)",
             "bench_impact/configs/impact_stanne_fidele_s15.cfg"),
            ("Percussion 2D (insert disque)", "configs/fdem_percussion.cfg"),
            ("Percussion 3D (insert sphère)",
             "configs/fdem3d_percussion.cfg"),
            ("Percussion 3D GBM Voronoï",
             "configs/fdem3d_voronoi_percussion.cfg"),
        ),
        groups=("Général", "Maillage", "Matériau", "Joints", "Contact",
                "Outil", "Corps et groupes", "Conditions aux limites",
                "Sorties"),
    ),
    Module(
        "labo", "Essais de laboratoire",
        "UCS, brésilien, triaxial — calibration Red Bohus.",
        templates=(
            ("UCS Bohus (platines)", "configs/cal_ucs_bohus.cfg"),
            ("Brésilien Bohus (disque)", "configs/cal_bts_bohus.cfg"),
            ("Triaxial Bohus GBM", "configs/triax_bohus_gbm.cfg"),
            ("Vérification traction FDEM",
             "configs/verify_fdem_tension.cfg"),
        ),
        groups=("Général", "Maillage", "Matériau", "Joints", "Contact",
                "Conditions aux limites", "Sorties"),
    ),
    Module(
        "expert", "Expert (tout)",
        "L'interface complète : tous les gabarits, tous les groupes.",
    ),
)


def by_key(key: str) -> Module:
    for m in MODULES:
        if m.key == key:
            return m
    return MODULES[-1]                     # expert par défaut


def sweep_cases(base_pairs: dict, sweep: Sweep, values: list[float]):
    """(suffixe, dict de clés modifiées) pour chaque valeur du balayage.
    PUR (aucun I/O) : testable sans interface."""
    out = []
    for v in values:
        tag = ("%g" % v).replace(".", "p")
        out.append((tag, dict(base_pairs, **sweep.make(v))))
    return out
