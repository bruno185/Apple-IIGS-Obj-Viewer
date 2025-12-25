# Build with Visual Studio (PowerShell)
$root = $PSScriptRoot
$build = Join-Path $root 'build_vs'
if (-not (Test-Path $build)) { New-Item -ItemType Directory -Path $build | Out-Null }
# Try several common VS generators (Visual Studio 2019/2022)
$gens = @("Visual Studio 17 2022","Visual Studio 16 2019")
$ok = $false
foreach ($g in $gens) {
    Write-Host "Trying generator: $g"
    cmake -S $root -B $build -G "$g" -A x64
    if ($LASTEXITCODE -eq 0) { $ok = $true; break }
}
if (-not $ok) { Write-Error "cmake generation failed (no supported VS generator found)"; exit 1 }
cmake --build $build --config Release
if ($LASTEXITCODE -ne 0) { Write-Error "build failed"; exit 1 }
Write-Host "Build complete: $build\bin\Release\viewer_win32.exe (or bin\Release)"