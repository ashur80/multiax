<#
.SYNOPSIS
  Build and test everything: CMake preset (C++ library, doctest, benchmarks), Python module, pytest.

.EXAMPLE
  powershell -File tools/check.ps1              # release preset
  powershell -File tools/check.ps1 -Preset debug
  powershell -File tools/check.ps1 -SkipPython  # C++ only
#>
param(
    [ValidateSet('release', 'debug')][string]$Preset = 'release',
    [switch]$SkipPython,
    [switch]$Bench
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

# llvm-mingw and Ninja from winget (in case this shell was started before the installation)
$pkgs = Join-Path $env:LOCALAPPDATA 'Microsoft\WinGet\Packages'
$llvm = Get-ChildItem -Path $pkgs -Filter 'MartinStorsjo.LLVM-MinGW.UCRT*' -Directory -ErrorAction SilentlyContinue |
        ForEach-Object { Get-ChildItem $_.FullName -Filter 'llvm-mingw-*' -Directory } | Select-Object -First 1
$ninja = Get-ChildItem -Path $pkgs -Filter 'Ninja-build.Ninja*' -Directory -ErrorAction SilentlyContinue | Select-Object -First 1
if ($llvm) { $env:PATH = (Join-Path $llvm.FullName 'bin') + ';' + $env:PATH }
if ($ninja) { $env:PATH = $ninja.FullName + ';' + $env:PATH }

function Invoke-Step([string]$name, [scriptblock]$cmd) {
    Write-Host "==> $name" -ForegroundColor Cyan
    & $cmd
    if ($LASTEXITCODE -ne 0) { throw "$name failed (exit $LASTEXITCODE)" }
}

Invoke-Step "configure ($Preset)" { cmake --preset $Preset }
Invoke-Step "build ($Preset)"     { cmake --build --preset $Preset }
Invoke-Step "C++ tests"           { ctest --preset $Preset }
if ($Bench) {
    Invoke-Step "benchmarks" { & "build/$Preset/bench/multiax_bench.exe" }
}
if (-not $SkipPython) {
    Invoke-Step "pip install -e ."  { python -m pip install --user --no-build-isolation -e . --quiet }
    Invoke-Step "Python tests"      { python -m pytest -q tests/python }
}
Write-Host "All checks passed." -ForegroundColor Green
