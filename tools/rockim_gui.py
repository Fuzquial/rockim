#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# rockim_gui.py — interface graphique de pilotage de rockim.
#
#   python tools/rockim_gui.py            (depuis la racine du projet)
#
# Trois zones :
#   * gauche  : les configs (configs/*.cfg), éditables et sauvegardables ;
#   * centre  : lancement (exe, dossier de sortie, threads OpenMP), journal
#               du run en direct, bouton "suite de vérification" ;
#   * droite  : résultats — courbe force-pénétration, historiques, coupe
#               médiane (fem3d/fdem 2D), pour n'importe quel dossier out_*.
#
# Dépendances : tkinter (standard) + numpy + matplotlib.
# ---------------------------------------------------------------------------
import os
import queue
import re
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

VERIF_CONFIGS = [
    ("FDEM 2D tension (grid)", "verify_fdem_tension.cfg"),
    ("FEM 2D onde de barre", "verify_fem_bar.cfg"),
    ("DEM 2D tension", "verify_dem_tension.cfg"),
    ("DEM3D tension", "verify_dem3d_tension.cfg"),
    ("FDEM3D tension (grid)", "verify_fdem3d_tension.cfg"),
    ("FDEM Voronoï tension", "verify_fdem_voronoi_tension.cfg"),
    ("FDEM3D Voronoï tension", "verify_fdem3d_voronoi_tension.cfg"),
    ("FEM3D compression DP", "verify_fem3d_dp.cfg"),
    ("FEM3D tension ft", "verify_fem3d_tension.cfg"),
    ("FEM3D Perzyna rate 1", "verify_fem3d_rate1.cfg"),
    ("FEM3D Perzyna rate 2", "verify_fem3d_rate2.cfg"),
]


def read_vtu(path, names):
    txt = open(path).read()
    m = re.search(r'<Points>.*?<DataArray[^>]*>\s*([^<]+)<', txt, re.S)
    pts = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    m = re.search(r'Name="connectivity"[^>]*>\s*([^<]+)<', txt)
    conn = np.fromstring(m.group(1), sep=" ").astype(int)
    m = re.search(r'Name="offsets"[^>]*>\s*([^<]+)<', txt)
    offs = np.fromstring(m.group(1), sep=" ").astype(int)
    ncell = offs[1] - offs[0] if len(offs) > 1 else len(conn)
    conn = conn.reshape(-1, ncell)
    data = {}
    for nm in names:
        m = re.search(r'Name="%s"[^>]*>\s*([^<]+)<' % nm, txt)
        if m:
            data[nm] = np.fromstring(m.group(1), sep=" ")
    return pts, conn, data


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("rockim — pilotage")
        self.geometry("1380x860")
        self.proc = None
        self.logq = queue.Queue()

        pw = ttk.PanedWindow(self, orient="horizontal")
        pw.pack(fill="both", expand=True)

        # ---- gauche : configs -------------------------------------------
        left = ttk.Frame(pw, width=430)
        pw.add(left, weight=1)
        ttk.Label(left, text="Configurations").pack(anchor="w", padx=4)
        self.cfgList = tk.Listbox(left, height=14, exportselection=False)
        self.cfgList.pack(fill="x", padx=4)
        self.cfgList.bind("<<ListboxSelect>>", self.load_cfg)
        self.cfgText = tk.Text(left, wrap="none", font=("Consolas", 9))
        self.cfgText.pack(fill="both", expand=True, padx=4, pady=4)
        fb = ttk.Frame(left)
        fb.pack(fill="x", padx=4, pady=2)
        ttk.Button(fb, text="Enregistrer", command=self.save_cfg).pack(
            side="left")
        ttk.Button(fb, text="Enregistrer sous…",
                   command=self.save_cfg_as).pack(side="left", padx=4)
        ttk.Button(fb, text="Rafraîchir la liste",
                   command=self.fill_cfgs).pack(side="right")

        # ---- centre : lancement + journal --------------------------------
        mid = ttk.Frame(pw, width=470)
        pw.add(mid, weight=1)
        row = ttk.Frame(mid)
        row.pack(fill="x", padx=4, pady=2)
        ttk.Label(row, text="Exécutable :").pack(side="left")
        self.exeVar = tk.StringVar(value=self.find_exe())
        ttk.Entry(row, textvariable=self.exeVar).pack(
            side="left", fill="x", expand=True, padx=4)
        ttk.Button(row, text="…", width=3, command=self.pick_exe).pack(
            side="left")
        row2 = ttk.Frame(mid)
        row2.pack(fill="x", padx=4, pady=2)
        ttk.Label(row2, text="Dossier de sortie :").pack(side="left")
        self.outVar = tk.StringVar(value="out_gui")
        ttk.Entry(row2, textvariable=self.outVar, width=18).pack(
            side="left", padx=4)
        ttk.Label(row2, text="Threads OMP :").pack(side="left", padx=(12, 0))
        self.thrVar = tk.StringVar(value="")
        ttk.Entry(row2, textvariable=self.thrVar, width=5).pack(side="left")
        row3 = ttk.Frame(mid)
        row3.pack(fill="x", padx=4, pady=4)
        self.runBtn = ttk.Button(row3, text="▶ Lancer la config",
                                 command=self.run_cfg)
        self.runBtn.pack(side="left")
        ttk.Button(row3, text="■ Arrêter", command=self.stop_run).pack(
            side="left", padx=6)
        ttk.Button(row3, text="✔ Suite de vérification",
                   command=self.run_verifs).pack(side="right")
        self.log = tk.Text(mid, wrap="word", font=("Consolas", 9),
                           state="disabled", bg="#111", fg="#ddd")
        self.log.pack(fill="both", expand=True, padx=4, pady=4)

        # ---- droite : résultats -----------------------------------------
        right = ttk.Frame(pw, width=480)
        pw.add(right, weight=2)
        rrow = ttk.Frame(right)
        rrow.pack(fill="x", padx=4, pady=2)
        ttk.Label(rrow, text="Dossier résultats :").pack(side="left")
        self.resVar = tk.StringVar(value="")
        self.resCombo = ttk.Combobox(rrow, textvariable=self.resVar,
                                     values=[], width=24)
        self.resCombo.pack(side="left", padx=4)
        ttk.Button(rrow, text="↻", width=3, command=self.fill_results).pack(
            side="left")
        rrow2 = ttk.Frame(right)
        rrow2.pack(fill="x", padx=4, pady=2)
        ttk.Button(rrow2, text="Courbe F–δ",
                   command=lambda: self.plot("fp")).pack(side="left")
        ttk.Button(rrow2, text="Historique",
                   command=lambda: self.plot("hist")).pack(side="left", padx=4)
        ttk.Button(rrow2, text="Coupe médiane",
                   command=lambda: self.plot("slice")).pack(side="left")
        self.fig = Figure(figsize=(6, 5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True,
                                         padx=4, pady=4)

        self.fill_cfgs()
        self.fill_results()
        self.after(100, self.poll_log)

    # ---------------- config management ----------------
    def find_exe(self):
        for c in ("rockim.exe", "build/rockim.exe", "build/Release/rockim.exe",
                  "rockim"):
            p = os.path.join(ROOT, c)
            if os.path.exists(p):
                return p
        return os.path.join(ROOT, "rockim.exe")

    def pick_exe(self):
        p = filedialog.askopenfilename(initialdir=ROOT)
        if p:
            self.exeVar.set(p)

    def fill_cfgs(self):
        self.cfgList.delete(0, "end")
        d = os.path.join(ROOT, "configs")
        for f in sorted(os.listdir(d)):
            if f.endswith(".cfg"):
                self.cfgList.insert("end", f)

    def cur_cfg(self):
        sel = self.cfgList.curselection()
        return self.cfgList.get(sel[0]) if sel else None

    def load_cfg(self, _=None):
        f = self.cur_cfg()
        if not f:
            return
        with open(os.path.join(ROOT, "configs", f), encoding="utf-8",
                  errors="replace") as fh:
            self.cfgText.delete("1.0", "end")
            self.cfgText.insert("1.0", fh.read())

    def save_cfg(self):
        f = self.cur_cfg()
        if not f:
            return
        with open(os.path.join(ROOT, "configs", f), "w",
                  encoding="ascii", errors="replace") as fh:
            fh.write(self.cfgText.get("1.0", "end-1c"))
        self.log_line(f"[gui] {f} enregistré")

    def save_cfg_as(self):
        p = filedialog.asksaveasfilename(
            initialdir=os.path.join(ROOT, "configs"),
            defaultextension=".cfg")
        if p:
            with open(p, "w", encoding="ascii", errors="replace") as fh:
                fh.write(self.cfgText.get("1.0", "end-1c"))
            self.fill_cfgs()

    # ---------------- running ----------------
    def log_line(self, s):
        self.log.configure(state="normal")
        self.log.insert("end", s.rstrip() + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def run_cfg(self):
        f = self.cur_cfg()
        if not f:
            messagebox.showinfo("rockim", "Sélectionne une config.")
            return
        self.save_cfg()
        cfg = os.path.join(ROOT, "configs", f)
        self.launch([self.exeVar.get(), cfg, self.outVar.get()])

    def launch(self, cmd, on_done=None):
        if self.proc is not None:
            messagebox.showinfo("rockim", "Un run est déjà en cours.")
            return
        env = os.environ.copy()
        if self.thrVar.get().strip():
            env["OMP_NUM_THREADS"] = self.thrVar.get().strip()
        self.log_line("[gui] $ " + " ".join(cmd))
        try:
            self.proc = subprocess.Popen(
                cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1)
        except OSError as e:
            self.log_line(f"[gui] échec du lancement : {e}")
            self.proc = None
            return
        self.runBtn.configure(state="disabled")

        def reader(p, q, done):
            for line in p.stdout:
                q.put(line)
            p.wait()
            q.put(f"[gui] terminé (code {p.returncode})\n")
            q.put(("__DONE__", done))

        threading.Thread(target=reader,
                         args=(self.proc, self.logq, on_done),
                         daemon=True).start()

    def poll_log(self):
        try:
            while True:
                item = self.logq.get_nowait()
                if isinstance(item, tuple) and item[0] == "__DONE__":
                    self.proc = None
                    self.runBtn.configure(state="normal")
                    self.fill_results()
                    if item[1]:
                        item[1]()
                else:
                    self.log_line(item)
        except queue.Empty:
            pass
        self.after(100, self.poll_log)

    def stop_run(self):
        if self.proc is not None:
            self.proc.terminate()
            self.log_line("[gui] arrêt demandé")

    def run_verifs(self):
        todo = list(VERIF_CONFIGS)
        results = []

        def next_one():
            if not todo:
                self.log_line("[gui] ---- SUITE DE VÉRIFICATION ----")
                for nm, verdict in results:
                    self.log_line(f"[gui]   {nm:28s} {verdict}")
                return
            nm, cfg = todo.pop(0)
            outd = "out_verif"
            logmark = len(self.log.get('1.0', 'end-1c'))

            def done():
                txt = self.log.get("1.0", "end-1c")[logmark:]
                verdict = ("PASS" if "PASS" in txt
                           else "FAIL" if "FAIL" in txt else "?")
                results.append((nm, verdict))
                next_one()

            self.launch([self.exeVar.get(),
                         os.path.join(ROOT, "configs", cfg), outd],
                        on_done=done)

        next_one()

    # ---------------- results ----------------
    def fill_results(self):
        outs = sorted(d for d in os.listdir(ROOT)
                      if d.startswith("out") and
                      os.path.isdir(os.path.join(ROOT, d)))
        self.resCombo.configure(values=outs)
        if outs and not self.resVar.get():
            self.resVar.set(outs[-1])

    def plot(self, kind):
        run = os.path.join(ROOT, self.resVar.get())
        if not os.path.isdir(run):
            messagebox.showinfo("rockim", "Choisis un dossier de résultats.")
            return
        self.fig.clf()
        try:
            if kind == "fp":
                self.plot_fp(run)
            elif kind == "hist":
                self.plot_hist(run)
            else:
                self.plot_slice(run)
        except Exception as e:
            messagebox.showerror("rockim", f"Tracé impossible : {e}")
            return
        self.canvas.draw()

    def hist(self, run):
        return np.genfromtxt(os.path.join(run, "history.csv"),
                             delimiter=",", names=True)

    def plot_fp(self, run):
        d = self.hist(run)
        ax = self.fig.add_subplot(111)
        if "toolFz" in d.dtype.names:                  # 3D
            pen = (d["toolZ"][0] - d["toolZ"]) * 1e3
            F = d["toolFz"] / 1e3
            ax.set_ylabel("F [kN]")
        elif "toolFy" in d.dtype.names:                # 2D
            pen = (d["toolY"][0] - d["toolY"]) * 1e3
            F = -d["toolFy"] / 1e3
            ax.set_ylabel("F [kN/m]")
        else:
            raise RuntimeError("pas de colonnes outil dans history.csv")
        ax.plot(pen, F, lw=1.2, color="#1e4b8c")
        ax.set_xlabel("déplacement de l'outil [mm]")
        ax.set_title("force–pénétration")
        ax.grid(alpha=0.3)

    def plot_hist(self, run):
        d = self.hist(run)
        names = [n for n in d.dtype.names if n != "t"]
        n = len(names)
        for k, nm in enumerate(names):
            ax = self.fig.add_subplot((n + 2) // 3, 3, k + 1)
            ax.plot(d["t"] * 1e3, d[nm], lw=0.9)
            ax.set_title(nm, fontsize=8)
            ax.tick_params(labelsize=7)
        self.fig.suptitle("history.csv (t en ms)", fontsize=9)
        self.fig.tight_layout()

    def plot_slice(self, run):
        files = sorted(f for f in os.listdir(run)
                       if re.match(r"(fem3d|fdem)_\d+\.vtu$", f))
        if not files:
            raise RuntimeError("aucun frame .vtu (fem3d/fdem)")
        pts, conn, d = read_vtu(os.path.join(run, files[-1]),
                                ["damage", "vonMises", "eroded"])
        ax = self.fig.add_subplot(111)
        champ = "damage" if "damage" in d and d["damage"].max() > 0 \
            else "vonMises"
        vals = d[champ]
        if conn.shape[1] == 4:                          # tets: tranche y~mi
            cz = pts[conn].mean(axis=1)
            ymid = 0.5 * (pts[:, 1].min() + pts[:, 1].max())
            dy = (pts[:, 1].max() - pts[:, 1].min()) / 20 + 1e-12
            sel = np.abs(cz[:, 1] - ymid) < dy
            from matplotlib.collections import PolyCollection
            polys = [pts[conn[k]][:, [0, 2]][[0, 1, 2]]
                     for k in np.where(sel)[0]
                     if not ("eroded" in d and d["eroded"][k] >= 1.0)]
            v = np.array([vals[k] for k in np.where(sel)[0]
                          if not ("eroded" in d and d["eroded"][k] >= 1.0)])
            pc = PolyCollection(polys, array=v, cmap="Reds"
                                if champ == "damage" else "viridis",
                                edgecolors="none")
            ax.add_collection(pc)
            self.fig.colorbar(pc, ax=ax, label=champ)
            ax.autoscale()
        else:                                           # triangles 2D
            from matplotlib.collections import PolyCollection
            polys = [pts[c][:, :2] for c in conn]
            pc = PolyCollection(polys, array=vals, cmap="Reds"
                                if champ == "damage" else "viridis",
                                edgecolors="none")
            ax.add_collection(pc)
            self.fig.colorbar(pc, ax=ax, label=champ)
            ax.autoscale()
        ax.set_aspect("equal")
        ax.set_title(f"{files[-1]} — {champ}")


if __name__ == "__main__":
    App().mainloop()
