#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_configs.py — genere les 9 configs des deux balayages de Wang et al.
# (2024) a partir de tunnel_ref_s5_lam1.cfg, et le script de lancement.
#
#   python tunnel_edz/tools/make_configs.py
#
# Balayage 1 (leur section 5.1) : in situ HYDROSTATIQUE 3, 4, 5, 6, 7 MPa.
# Balayage 2 (leur section 5.2) : sigma_h = 5 MPa FIXE, sigma_v = 5/lambda
#   pour lambda = 0,5 / 0,75 / 1,0 / 1,25 / 1,5 (soit 10 / 6,67 / 5 / 4 / 3,33).
# Le cas (5 MPa, lambda = 1) est commun aux deux : il n'est ecrit qu'une fois.
#
# Sortie : tunnel_edz/configs/sweep/*.cfg + run_sweep.cmd
# ---------------------------------------------------------------------------
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CFGDIR = os.path.join(HERE, "..", "configs")
OUTDIR = os.path.join(CFGDIR, "sweep")
REF = os.path.join(CFGDIR, "tunnel_ref_s5_lam1.cfg")

# (nom, sigma_h [MPa], sigma_v [MPa], commentaire)
CASES = []
for s0 in (3.0, 4.0, 6.0, 7.0):                       # 5 = le cas de reference
    CASES.append((f"tunnel_s{int(s0)}_lam1", s0, s0,
                  f"balayage in situ (section 5.1) : sigma0 = {s0:g} MPa "
                  f"hydrostatique"))
for lam, sv in ((0.5, 10.0), (0.75, 6.67), (1.25, 4.0), (1.5, 3.33)):
    tag = f"{lam:g}".replace(".", "p")
    CASES.append((f"tunnel_lam{tag}", 5.0, sv,
                  f"balayage lambda (section 5.2) : lambda = {lam:g}, "
                  f"sigma_h = 5 MPa, sigma_v = {sv:g} MPa"))


def main():
    with open(REF) as f:
        ref = f.read().splitlines()
    # on jette l'en-tete de la reference (elle parle du cas 5 MPa) et on en
    # ecrit une propre par cas ; le corps garde ses commentaires en ligne.
    i = 0
    while i < len(ref) and (ref[i].startswith("#") or not ref[i].strip()):
        i += 1
    ref = ref[i:]
    os.makedirs(OUTDIR, exist_ok=True)
    runs = ["@echo off",
            "rem genere par tunnel_edz/tools/make_configs.py — a lancer depuis",
            "rem la RACINE du depot, un run a la fois (machine idle).",
            "rem Le cas de reference (5 MPa, lambda = 1) est a part :",
            "rem   rockim_tun.exe tunnel_edz\\configs\\tunnel_ref_s5_lam1.cfg out_tun_s5"]
    for name, sh, sv, note in CASES:
        out = ["# " + "-" * 73,
               "# " + note,
               "# Genere par tunnel_edz/tools/make_configs.py depuis",
               "# tunnel_edz/configs/tunnel_ref_s5_lam1.cfg — ne pas editer a",
               "# la main : reprendre la reference puis relancer le generateur.",
               "# EXIGE les patchs 1 (in situ) et 2 (excavation).",
               "#",
               f"#   rockim_tun.exe tunnel_edz\\configs\\sweep\\{name}.cfg "
               f"out_{name}",
               "#   (depuis la RACINE du depot : meshFile est relatif au cwd)",
               "# " + "-" * 73]
        for line in ref:
            k = line.split("=")[0].strip()
            if k == "insituSh":
                out.append(f"insituSh = {sh * 1e6:g}         # {sh:g} MPa")
            elif k == "insituSv":
                out.append(f"insituSv = {sv * 1e6:g}         # {sv:g} MPa "
                           f"(lambda = {sh / sv:.3g})")
            elif k == "frames":
                out.append("frames = 6                  # balayage : 6 trames "
                           "suffisent")
            else:
                out.append(line)
        path = os.path.join(OUTDIR, name + ".cfg")
        with open(path, "w") as f:
            f.write("\n".join(out) + "\n")
        print("ecrit :", os.path.relpath(path, os.path.join(HERE, "..", "..")))
        runs.append(f"rockim_tun.exe tunnel_edz\\configs\\sweep\\{name}.cfg "
                    f"out_{name}")
    runs.append("rem depouillement :")
    runs.append("rem   for %%d in (out_tun_* out_tunnel_*) do "
                "python tunnel_edz\\tools\\edz_metrics.py %%d")
    with open(os.path.join(OUTDIR, "run_sweep.cmd"), "w") as f:
        f.write("\r\n".join(runs) + "\r\n")
    print(f"{len(CASES)} configs + run_sweep.cmd")


if __name__ == "__main__":
    main()
