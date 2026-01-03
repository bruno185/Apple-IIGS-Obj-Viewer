# Cleanup script to remove session-context files and caches, then commit
$toRemove = @(
  'session_context.md',
  'scripts/update_session_context.py',
  'scripts/install_git_hook.py',
  'scripts/validate_scripts.py',
  'scripts/hooks/post_commit.py',
  'scripts/hooks/post-commit.wrapper.sh',
  'scripts/hooks/post-commit.example',
  'scripts/update_session_context.ps1',
  'scripts/install_git_hook.ps1',
  'scripts/validate_scripts.ps1'
)

foreach ($p in $toRemove) {
    if (Test-Path $p) {
        try { git rm -f -- $p; Write-Host "Removed: $p" }
        catch { Write-Host "git rm failed for: $p" }
    } else {
        Write-Host "Not found: $p"
    }
}

# Remove any __pycache__ under scripts
$pycaches = Get-ChildItem -Path scripts -Directory -Recurse -Force -ErrorAction SilentlyContinue | Where-Object { $_.Name -eq '__pycache__' }
foreach ($d in $pycaches) {
    try { git rm -r -f -- $d.FullName; Write-Host "Removed cache: $($d.FullName)" }
    catch { Write-Host "git rm failed for cache: $($d.FullName)" }
}

# Stage all other changes and commit if something staged
git add -A
$staged = (git diff --cached --name-only)
if ($staged) {
    git commit -m 'Remove session context scripts and artifacts (user requested)'
    Write-Host 'Committed deletions.'
} else {
    Write-Host 'No changes to commit.'
}
