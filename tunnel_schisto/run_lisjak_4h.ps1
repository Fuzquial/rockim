# ---------------------------------------------------------------------------
# run_lisjak_4h.ps1 — balayage du pendage, methode de Lisjak, reglage 4 h 30.
# A lancer DANS VOTRE terminal, depuis la racine de rockim_f2 :
#
#   powershell -ExecutionPolicy Bypass -File tunnel_schisto\run_lisjak_4h.ps1
#
# Trois runs SEQUENTIELS (45, puis 0, puis 90 deg), 14 threads, un job a la
# fois : refuse de demarrer si un rockim tourne deja. Journal par run dans
# out_lisjakXX_4h\run.log ; history.csv y est vide toutes les ~2000 lignes,
# donc lisible en cours de run. Arret propre : Ctrl+C dans ce terminal.
#
# Cout mesure au dry-run (2026-09-02) : dt = 3,22e-6 s, 217 730 pas,
# ~4 h 30 par pendage sur 14 threads, ~13 h 30 au total.
# ---------------------------------------------------------------------------
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path))

$busy = Get-Process | Where-Object { $_.ProcessName -match '^rockim' }
if ($busy) {
    Write-Host "Un rockim tourne deja (PID $($busy.Id -join ', ')) : un job a la fois. Abandon." -ForegroundColor Red
    exit 1
}
if (-not (Test-Path .\rockim_f2i.exe)) { Write-Host "rockim_f2i.exe introuvable" -ForegroundColor Red; exit 1 }

$env:OMP_NUM_THREADS = "14"
foreach ($a in "45", "00", "90") {
    $deck = "tunnel_schisto\tunnel_lisjak${a}_4h.cfg"
    $out  = "out_lisjak${a}_4h"
    $mesh = "meshes\tunnel_hs_bed${a}_t06.msh"
    if (-not (Test-Path $deck)) { Write-Host "deck manquant : $deck" -ForegroundColor Red; exit 1 }
    if (-not (Test-Path $mesh)) { Write-Host "maillage manquant : $mesh" -ForegroundColor Red; exit 1 }
    New-Item -ItemType Directory -Force $out | Out-Null
    $t0 = Get-Date
    Write-Host ("[{0:HH:mm}] pendage {1} deg -> {2}" -f $t0, [int]$a, $out) -ForegroundColor Cyan
    & .\rockim_f2i.exe $deck $out 2>&1 | Tee-Object -FilePath "$out\run.log"
    $dtm = (Get-Date) - $t0
    Write-Host ("[{0:HH:mm}] termine en {1:h\ h\ mm\ min}, code {2}" -f (Get-Date), $dtm, $LASTEXITCODE) -ForegroundColor Cyan
    if ($LASTEXITCODE -ne 0) { Write-Host "run en erreur : arret du balayage" -ForegroundColor Red; exit $LASTEXITCODE }
}
Write-Host "Balayage termine. Depouillement :" -ForegroundColor Green
Write-Host "  python ..\rockim\rockim_p1\tunnel_edz\tools\edz_metrics.py       out_lisjak45_4h"
Write-Host "  python ..\rockim\rockim_p1\tunnel_edz\tools\crack_orientation.py out_lisjak45_4h"
