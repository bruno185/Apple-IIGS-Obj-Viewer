#!/usr/bin/env python3
"""Install the example git hook wrapper into .git/hooks/post-commit
Usage: python scripts/install_git_hook.py
"""
from pathlib import Path
import shutil
import stat

ROOT = Path(__file__).resolve().parent.parent
HOOK_SRC = ROOT / 'scripts' / 'hooks' / 'post-commit.wrapper.sh'
HOOK_DST = ROOT / '.git' / 'hooks' / 'post-commit'


def install():
    if not (ROOT / '.git').exists():
        print('.git directory not found. Are you in the repo root?')
        return 1
    dest = HOOK_DST
    shutil.copy2(HOOK_SRC, dest)
    # make executable
    try:
        st = dest.stat()
        dest.chmod(st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except Exception:
        pass
    print(f'Installed post-commit hook to {dest}')
    print('Note: The hook will append the last commit info to session_context.md and commit it. Ensure this workflow is acceptable for your CI.')
    return 0


if __name__ == '__main__':
    raise SystemExit(install())
