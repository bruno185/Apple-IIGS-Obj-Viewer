#!/usr/bin/env python3
"""Hook helper: append last commit info to session_context.md and commit it safely.
Intended to be called by a shell wrapper placed in .git/hooks/post-commit
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent.parent
SESSION_FILE = ROOT / 'session_context.md'


def last_commit_info():
    # Get last commit short hash, author, subject
    p = subprocess.run(['git', 'log', '-1', '--pretty=format:%h | %an | %s'], cwd=str(ROOT), capture_output=True, text=True)
    if p.returncode != 0:
        return None
    return p.stdout.strip()


def last_commit_message_contains_marker():
    p = subprocess.run(['git', 'log', '-1', '--pretty=%B'], cwd=str(ROOT), capture_output=True, text=True)
    if p.returncode != 0:
        return False
    return '[session context]' in p.stdout


def append_and_commit(info: str):
    from datetime import datetime, timezone
    date = datetime.now(timezone.utc).isoformat(timespec='seconds')
    line = f"{date} | {info}"
    try:
        with open(SESSION_FILE, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')
    except Exception as e:
        print('Failed to append to session_context.md:', e, file=sys.stderr)
        return False
    # Commit with a marker message so the hook won't recurse
    commit_msg = f"[session context] auto-update: {info}"
    p = subprocess.run(['git', 'add', str(SESSION_FILE)], cwd=str(ROOT))
    if p.returncode != 0:
        print('git add failed', file=sys.stderr)
        return False
    p = subprocess.run(['git', 'commit', '-m', commit_msg], cwd=str(ROOT))
    if p.returncode != 0:
        print('git commit failed', file=sys.stderr)
        return False
    return True


def main():
    if last_commit_message_contains_marker():
        # avoid recursion
        sys.exit(0)
    info = last_commit_info()
    if not info:
        sys.exit(0)
    ok = append_and_commit(info)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
