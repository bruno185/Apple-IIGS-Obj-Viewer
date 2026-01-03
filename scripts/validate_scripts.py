#!/usr/bin/env python3
"""Quick validator for Python scripts in scripts/ - checks syntax by compiling. Returns exit code 0 on success."""
import glob
import py_compile
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PATTERN = str(ROOT / 'scripts' / '*.py')

failed = False
for path in glob.glob(PATTERN):
    try:
        py_compile.compile(path, doraise=True)
        print(f'OK: {path}')
    except py_compile.PyCompileError as e:
        print(f'ERROR in {path}: {e.msg}')
        failed = True

if failed:
    sys.exit(1)
else:
    sys.exit(0)
