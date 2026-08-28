# ---------------------------------------------------------------------------
# pause_run.ps1 - SUSPENDRE / REPRENDRE un run rockim en cours, sans le perdre.
#
#   powershell -File tools\pause_run.ps1 -ProcId 33336 -Action pause
#   powershell -File tools\pause_run.ps1 -ProcId 33336 -Action reprendre
#   powershell -File tools\pause_run.ps1 -ProcId 33336 -Action etat
#
# POURQUOI CE SCRIPT. rockim n'a AUCUN mecanisme de reprise sur fichier
# (verifie le 2026-08-28 : aucun restart / checkpoint dans src/ ni dans la
# doc). Tuer un run, c'est donc perdre tout ce qui n'a pas ete ecrit en frame,
# et le relancer repart de t = 0. La seule facon de mettre en pause est donc
# de GELER LE PROCESSUS : NtSuspendProcess fige les threads en place, la
# memoire reste allouee, et NtResumeProcess repart au pas exact ou l'on s'est
# arrete. Mesure de controle : 0,00 s de CPU consomme pendant la suspension.
#
# CE QUE CA NE PROTEGE PAS : un processus suspendu meurt avec la session
# Windows (deconnexion, redemarrage, coupure de courant). La pause sert a
# liberer la machine quelques heures, PAS a traverser un reboot. Pour cela il
# faudrait un vrai checkpoint dans le solveur - chantier ouvert.
#
# La RAM reste occupee pendant la pause (elle peut partir en swap : la reprise
# est alors lente les premieres secondes, c'est normal).
#
# NOTE : le parametre s'appelle ProcId et non Pid, $PID etant une variable
# automatique reservee de PowerShell. Et ce fichier est en ASCII PUR : les
# accents et tirets longs en .ps1 arrivent en mojibake sur cette machine et
# cassent le parseur (piege documente du projet).
# ---------------------------------------------------------------------------
param(
  [Parameter(Mandatory=$true)][int]$ProcId,
  [ValidateSet("pause","reprendre","etat")][string]$Action = "etat"
)

$src = @"
using System;
using System.Runtime.InteropServices;
public static class RockimSusp {
  [DllImport("ntdll.dll")] public static extern uint NtSuspendProcess(IntPtr h);
  [DllImport("ntdll.dll")] public static extern uint NtResumeProcess(IntPtr h);
}
"@
if (-not ("RockimSusp" -as [type])) { Add-Type -TypeDefinition $src }

try { $p = Get-Process -Id $ProcId -ErrorAction Stop }
catch { Write-Host "PID $ProcId introuvable : le run n'est plus la." ; exit 1 }

if ($p.ProcessName -notlike "rockim*") {
  Write-Host "PID $ProcId = '$($p.ProcessName)' : ce n'est pas un run rockim. Arret."
  exit 1
}

function Mesure-CPU {
  $a = (Get-Process -Id $ProcId).TotalProcessorTime
  Start-Sleep -Seconds 3
  $b = (Get-Process -Id $ProcId).TotalProcessorTime
  return ($b - $a).TotalSeconds
}

if ($Action -eq "pause") {
  [void][RockimSusp]::NtSuspendProcess($p.Handle)
  $c = Mesure-CPU
  if ($c -lt 0.05) {
    Write-Host "SUSPENDU (PID $ProcId) : 0 CPU, la machine est libre."
  } else {
    Write-Host ("ECHEC : le processus consomme encore {0:N2} s de CPU." -f $c)
  }
}
elseif ($Action -eq "reprendre") {
  [void][RockimSusp]::NtResumeProcess($p.Handle)
  $c = Mesure-CPU
  if ($c -gt 1.0) {
    Write-Host ("REPRIS (PID $ProcId) : {0:N1} coeurs actifs." -f ($c/3))
  } else {
    Write-Host "ATTENTION : toujours gele. Relancer 'reprendre' (les suspensions s'empilent)."
  }
}
else {
  $c = Mesure-CPU
  if ($c -lt 0.05) { $etat = "SUSPENDU" }
  else { $etat = "EN COURS ({0:N1} coeurs)" -f ($c/3) }
  Write-Host ("PID $ProcId : $etat, RAM {0:N0} Mo" -f ($p.WorkingSet64/1MB))
}
