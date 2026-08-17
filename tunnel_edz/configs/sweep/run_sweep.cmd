@echo off
rem genere par tunnel_edz/tools/make_configs.py — a lancer depuis
rem la RACINE du depot, un run a la fois (machine idle).
rem Le cas de reference (5 MPa, lambda = 1) est a part :
rem   rockim_tun.exe tunnel_edz\configs\tunnel_ref_s5_lam1.cfg out_tun_s5
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_s3_lam1.cfg out_tunnel_s3_lam1
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_s4_lam1.cfg out_tunnel_s4_lam1
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_s6_lam1.cfg out_tunnel_s6_lam1
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_s7_lam1.cfg out_tunnel_s7_lam1
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_lam0p5.cfg out_tunnel_lam0p5
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_lam0p75.cfg out_tunnel_lam0p75
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_lam1p25.cfg out_tunnel_lam1p25
rockim_tun.exe tunnel_edz\configs\sweep\tunnel_lam1p5.cfg out_tunnel_lam1p5
rem depouillement :
rem   for %%d in (out_tun_* out_tunnel_*) do python tunnel_edz\tools\edz_metrics.py %%d
