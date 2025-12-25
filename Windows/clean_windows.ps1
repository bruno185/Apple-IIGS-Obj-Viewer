Param(
    [switch]$WhatIf,
    [switch]$Force
)

<#
Safe cleanup script for Windows/build_vs generated files.
- Shows the list of files/folders to remove (if any)
- Asks for confirmation (unless -Force is used)
- Moves files/folders to the Recycle Bin using Microsoft.VisualBasic.FileIO
Examples:
  .\clean_windows.ps1        # interactive confirmation
  .\clean_windows.ps1 -WhatIf # show what would be removed
  .\clean_windows.ps1 -Force  # non-interactive, perform deletion
#>

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$targets = @()

function Add-IfExists($relative) {
    $full = Join-Path $root $relative
    if (Test-Path $full) { $script:targets += $full }
}

# Common generated files and directories to remove
$toCheck = @(
    'build_vs\GS3DpViewer.sln',
    'build_vs\ALL_BUILD.vcxproj',
    'build_vs\ALL_BUILD.vcxproj.filters',
    'build_vs\ZERO_CHECK.vcxproj',
    'build_vs\ZERO_CHECK.vcxproj.filters',
    'build_vs\viewer.vcxproj',
    'build_vs\viewer.vcxproj.filters',
    'build_vs\viewer_win32.vcxproj',
    'build_vs\viewer_win32.vcxproj.filters',
    'build_vs\CMakeCache.txt',
    'build_vs\cmake_install.cmake',
    'build_vs\CMakeFiles',
    'build_vs\bin',
    'build_vs\viewer.dir',
    'build_vs\viewer_win32.dir',
    'build_vs\ALL_BUILD.dir',
    'build_vs\x64'
)

foreach ($r in $toCheck) { Add-IfExists $r }

# Also add common temp files under build_vs (tlog, recipe, lastbuildstate, exe, obj)
$extras = Get-ChildItem -Path (Join-Path $root 'build_vs') -Recurse -Force -ErrorAction SilentlyContinue |
    Where-Object { $_.Extension -in '.tlog','.recipe','.lastbuildstate','.exe','.obj' } | ForEach-Object { $_.FullName }

$targets += $extras

# Deduplicate
$targets = $targets | Sort-Object -Unique

if ($targets.Count -eq 0) {
    Write-Host "Aucun fichier généré trouvé dans build_vs. Rien à faire." -ForegroundColor Green
    exit 0
}

Write-Host "Les fichiers/dossiers suivants seront envoyés à la Corbeille:" -ForegroundColor Yellow
$targets | ForEach-Object { Write-Host "  $_" }

if ($WhatIf) { Write-Host "Mode simulation : aucun changement effectué" -ForegroundColor Cyan; exit 0 }

if (-not $Force) {
    $answer = Read-Host "Confirmez-vous la suppression et l'envoi à la Corbeille ? (O/N)"
    if ($answer -notin 'O','o','Y','y') { Write-Host "Annulé par l'utilisateur." -ForegroundColor Red; exit 0 }
}

# Use Microsoft.VisualBasic.FileIO to send items to Recycle Bin
Add-Type -AssemblyName Microsoft.VisualBasic

foreach ($p in $targets) {
    try {
        if (Test-Path $p -PathType Container) {
            [Microsoft.VisualBasic.FileIO.FileSystem]::DeleteDirectory($p, [Microsoft.VisualBasic.FileIO.UIOption]::OnlyErrorDialogs, [Microsoft.VisualBasic.FileIO.RecycleOption]::SendToRecycleBin)
            Write-Host "Déplacé vers Corbeille : $p" -ForegroundColor Green
        } else {
            [Microsoft.VisualBasic.FileIO.FileSystem]::DeleteFile($p, [Microsoft.VisualBasic.FileIO.UIOption]::OnlyErrorDialogs, [Microsoft.VisualBasic.FileIO.RecycleOption]::SendToRecycleBin)
            Write-Host "Déplacé vers Corbeille : $p" -ForegroundColor Green
        }
    } catch {
        Write-Warning "Échec pour $p : $_"
    }
}

Write-Host "Nettoyage terminé." -ForegroundColor Green
