# ---------------------------------------------------------------------------
# run_solutions.ps1 — enchaine des decks de la suite « traversee du litage ».
# A lancer DANS VOTRE terminal, depuis la racine de rockim_f2 :
#
#   powershell -ExecutionPolicy Bypass -File tunnel_schisto\run_solutions.ps1 S4_ratios1 S1_G030 S3_lambda130
#
# Un job a la fois (refuse si un rockim tourne), 14 threads, journal dans
# out_<deck>\run.log, history.csv lisible en cours de run. Ctrl+C = arret propre.
# Cout mesure : ~4 h 11 par deck (rockim_f2j.exe, 117 132 triangles).
# ---------------------------------------------------------------------------
param([Parameter(ValueFromRemainingArguments = $true)][string[]]$decks)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path))
if (-not $decks) { Write-Host "usage : run_solutions.ps1 <deck> [<deck> ...]   (noms sans .cfg, dans tunnel_schisto\)"; exit 1 }
$busy = Get-Process | Where-Object { $_.ProcessName -match '^rockim' }
if ($busy) { Write-Host "Un rockim tourne deja (PID $($busy.Id -join ', ')) : un job a la fois. Abandon." -ForegroundColor Red; exit 1 }
if (-not (Test-Path .\rockim_f2j.exe)) { Write-Host "rockim_f2j.exe introuvable" -ForegroundColor Red; exit 1 }
$env:OMP_NUM_THREADS = "14"
foreach ($d in $decks) {
    $deck = "tunnel_schisto\$d.cfg"; $out = "out_$d"
    if (-not (Test-Path $deck)) { Write-Host "deck manquant : $deck" -ForegroundColor Red; exit 1 }
    New-Item -ItemType Directory -Force $out | Out-Null
    $t0 = Get-Date
    Write-Host ("[{0:HH:mm}] {1} -> {2}" -f $t0, $d, $out) -ForegroundColor Cyan
    & .\rockim_f2j.exe $deck $out 2>&1 | Tee-Object -FilePath "$out\run.log"
    Write-Host ("[{0:HH:mm}] termine en {1:h\ h\ mm\ min}, code {2}" -f (Get-Date), ((Get-Date) - $t0), $LASTEXITCODE) -ForegroundColor Cyan
    if ($LASTEXITCODE -ne 0) { Write-Host "run en erreur : arret" -ForegroundColor Red; exit $LASTEXITCODE }
    Write-Host "  depouillement : python tunnel_schisto\tools\edz_sectors.py $out --dip 45 ; python tunnel_schisto\tools\joint_state_stats.py $out --dip 45"
}
Write-Host "Termine." -ForegroundColor Green
