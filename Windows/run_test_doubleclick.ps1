$exe = Join-Path $PSScriptRoot "build_vs\bin\Release\viewer_win32.exe"
Write-Host "Starting: $exe"
$p = Start-Process -FilePath $exe -PassThru -ErrorAction SilentlyContinue
Start-Sleep -Seconds 3
if ($p -and (Get-Process -Id $p.Id -ErrorAction SilentlyContinue)) {
    Write-Host "Process still running (Id=$($p.Id))"
    Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
} else {
    Write-Host "Process exited during startup"
    Write-Host "Last log lines:"
    Get-Content "$env:TEMP\viewer_win32.log" -Tail 40 | ForEach-Object { Write-Host $_ }
}
Write-Host "Done"