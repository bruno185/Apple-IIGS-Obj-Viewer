#!/bin/sh
# Wrapper to execute the Python hook helper from Git hooks
# Copy this file to .git/hooks/post-commit (or let the installer do it)
PYTHON=python
if command -v python3 >/dev/null 2>&1; then
  PYTHON=python3
fi
exec "$PYTHON" "$(git rev-parse --show-toplevel)/scripts/hooks/post_commit.py"
