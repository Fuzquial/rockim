# ---------------------------------------------------------------------------
# build.ps1 — construction de rockim_g0 par CMake + Ninja (mesure D1a).
# Remplace build_f2.cmd, dont les drapeaux sont reproduits a l'identique dans
# CMakeLists.txt (voir l'en-tete du CMakeLists : /nologo /std:c++17 /O2 /EHsc
# /openmp /D_USE_MATH_DEFINES /DNOMINMAX, CRT statique /MT).
#
#   powershell -ExecutionPolicy Bypass -File tools\build.ps1              # Release
#   powershell -ExecutionPolicy Bypass -File tools\build.ps1 -Clean       # reconfigure
#
# cmake.exe et ninja.exe sont ceux livres par Visual Studio 2022 (aucune
# installation supplementaire). Le resultat est build\rockim.exe.
# ---------------------------------------------------------------------------
param(
  [string]$BuildType = "Release",
  [string]$BuildDir  = "build",
  [int]$Jobs         = 4,
  [switch]$Clean
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot

$vsRoot = "C:\Program Files\Microsoft Visual Studio\2022\Community"
$cmake  = Join-Path $vsRoot "Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
$ninja  = Join-Path $vsRoot "Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe"
$vcvars = Join-Path $vsRoot "VC\Auxiliary\Build\vcvars64.bat"
foreach ($p in @($cmake, $ninja, $vcvars)) {
  if (-not (Test-Path $p)) { throw "introuvable : $p" }
}

$bd = Join-Path $root $BuildDir
if ($Clean -and (Test-Path $bd)) { Remove-Item -Recurse -Force $bd }

# vcvars64 exporte l'environnement du compilateur ; on l'importe dans CETTE
# session PowerShell pour que cmake trouve cl.exe.
$lines = & cmd /c "`"$vcvars`" >nul && set"
foreach ($l in $lines) {
  if ($l -match '^([^=]+)=(.*)$') { Set-Item -Path "Env:$($Matches[1])" -Value $Matches[2] }
}

# Arguments passes en TABLEAU : en continuation de ligne (backtick), PowerShell
# 5.1 ne developpe pas toujours "$var" dans un bareword commencant par '-'
# (constate le 2026-09-05 : CMAKE_BUILD_TYPE valait litteralement '$BuildType',
# et ninja refusait le fichier de regles genere).
$cfgArgs = @(
  "-S", $root, "-B", $bd, "-G", "Ninja",
  "-DCMAKE_BUILD_TYPE=$BuildType",
  "-DCMAKE_MAKE_PROGRAM=$ninja",
  "-DCMAKE_CXX_COMPILER=cl"
)
& $cmake @cfgArgs
if ($LASTEXITCODE -ne 0) { throw "cmake configure : code $LASTEXITCODE" }

& $cmake --build $bd --parallel $Jobs
if ($LASTEXITCODE -ne 0) { throw "cmake build : code $LASTEXITCODE" }

$exe = Join-Path $bd "rockim.exe"
if (-not (Test-Path $exe)) { throw "rockim.exe absent apres build" }
$h = (Get-FileHash $exe -Algorithm SHA256).Hash.ToLower()
Write-Output "rockim.exe : $exe"
Write-Output "sha256     : $h"
Write-Output "taille     : $((Get-Item $exe).Length) octets"
