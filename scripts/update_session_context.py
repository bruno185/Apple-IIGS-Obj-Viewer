#!/usr/bin/env python3
"""Append a one-line session summary to session_context.md and optionally commit it.
Usage: python scripts/update_session_context.py --message "Short summary" --files "file1,file2" --status WIP
"""
import argparse
from datetime import datetime, timezone
import getpass
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SESSION_FILE = ROOT / 'session_context.md'


def run_git_commit(message: str) -> bool:
    try:
        subprocess.run(['git', 'add', str(SESSION_FILE)], check=True, cwd=str(ROOT))
        subprocess.run(['git', 'commit', '-m', f"[session context] {message}"], check=True, cwd=str(ROOT))
        print('Session context appended and committed.')
        return True
    except subprocess.CalledProcessError:
        print('Appended to session_context.md but git commit failed (git not available or not a repo).', file=sys.stderr)
        return False


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-m', '--message', help='Brief session summary (one short sentence)')
    p.add_argument('-f', '--files', help='Comma-separated relevant files', default='')
    p.add_argument('-s', '--status', help='Optional status (WIP, done, review)', default='WIP')
    args = p.parse_args()

    message = args.message
    if not message:
        try:
            message = input('Brief session summary (one short sentence): ').strip()
        except (KeyboardInterrupt, EOFError):
            print('\nAborted.')
            sys.exit(1)

    if not message:
        print('No message provided; aborting.')
        sys.exit(1)

    date = datetime.now(timezone.utc).isoformat(timespec='seconds')
    author = getpass.getuser()
    files = args.files or ''
    status = args.status or 'WIP'

    line = f"{date} | {author} | {message} | {status} | {files} |"
    try:
        with open(SESSION_FILE, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')
        print(f"Appended to {SESSION_FILE}")
    except Exception as e:
        print('Failed to append to session_context.md:', e, file=sys.stderr)
        sys.exit(1)

    # Try to commit, but do not fail if git is not available
    run_git_commit(message)


if __name__ == '__main__':
    main()
