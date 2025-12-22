#!/usr/bin/env python3
"""
Generate a markdown block of important functions with file:line references
and insert it into README.md between markers <!-- FUNC_LIST_START --> and <!-- FUNC_LIST_END -->.

Usage:
  python scripts/generate_function_refs.py GS3Dp.cc README.md

The script looks for each function name in the source file and records the first line
containing the name and parentheses. It writes a list like:
- 🔧 `signature` — `GS3Dp.cc:LINE`

This is a pragmatic tool and assumes prototypes/definitions appear on single lines.
"""
import re
import sys
from pathlib import Path

FUNC_NAMES = [
    'painter_newell_sancha_fast', 'painter_newell_sancha', 'dumpFaceEquationsCSV', 'createModel3D', 'destroyModel3D', 'loadModel3D',
    'computeModelBoundingSphere', 'computeDistanceFromBoundingSphere', 'getObserverParams',
    'processModelFast', 'processModelWireframe', 'readVertices', 'readFaces_model',
    'projectTo2D', 'calculateFaceDepths', 'computeDistanceToFit', 'autoScaleModel',
    'revertAutoScaleModel', 'backupModelCoords', 'freeBackupModelCoords', 'fitModelToView',
    'drawPolygons', 'main'
]

SIG_RE = re.compile(r"^\s*([\w\*\s]+?)\s+([A-Za-z_][A-Za-z0-9_]*)\s*\((.*?)\)\s*([;{]?)\s*$")


def find_signature(lines, name):
    # First try: match a clean prototype/definition line where the function name is the identifier
    for i, line in enumerate(lines):
        m = SIG_RE.match(line.strip())
        if m and m.group(2) == name:
            # Ensure this is a real prototype/definition (return type has letters), not a bullet/comment
            ret = m.group(1).strip()
            if not re.search(r"[A-Za-z]", ret):
                # skip matches where "return type" is not a true type (e.g., bullet lines like '* loadModel3D(...)')
                continue
            func = m.group(2).strip()
            args = m.group(3).strip()
            signature = f"{ret} {func}({args})"
            return i + 1, signature

    # Second try: find any line that contains the name followed by '(' and try to expand across a small window
    for i, line in enumerate(lines):
        if re.search(r"\b" + re.escape(name) + r"\s*\(", line):
            start = max(0, i - 3)
            end = min(len(lines) - 1, i + 5)
            window = ' '.join(l.strip() for l in lines[start:end + 1])
            # try to extract a proper signature from the window
            m = SIG_RE.search(window)
            if m and m.group(2) == name:
                ret = m.group(1).strip()
                func = m.group(2).strip()
                args = m.group(3).strip()
                signature = f"{ret} {func}({args})"
                return i + 1, signature

            # fallback: extract substring from first occurrence of name up to the closing ')' and optional '{' or ';'
            m2 = re.search(re.escape(name) + r"\s*\([^\)]*\)\s*[;{]?", window)
            if m2:
                signature = m2.group(0).strip().rstrip('{').rstrip(';')
                return i + 1, signature

            # last resort: return the trimmed line
            signature = line.strip().rstrip('{').rstrip(';')
            return i + 1, signature

    return None, None


def generate_block(src_path):
    lines = src_path.read_text(encoding='utf-8').splitlines()
    md_lines = []
    for name in FUNC_NAMES:
        ln, sig = find_signature(lines, name)
        if ln and sig:
            md_lines.append(f"- 🔧 `{sig}` — `{src_path.name}:{ln}`")
        else:
            md_lines.append(f"- 🔧 `{name}(...)` — `MISSING`")
    return "\n".join(md_lines)


def replace_in_readme(readme_path, block_text):
    text = readme_path.read_text(encoding='utf-8')
    start = '<!-- FUNC_LIST_START -->'
    end = '<!-- FUNC_LIST_END -->'
    if start in text and end in text:
        pre, rest = text.split(start, 1)
        _, post = rest.split(end, 1)
        new_text = pre + start + '\n' + block_text + '\n' + end + post
        readme_path.write_text(new_text, encoding='utf-8')
        return True
    else:
        print("Markers not found in README; no changes made.")
        return False


def main():
    if len(sys.argv) < 3:
        print("Usage: generate_function_refs.py <source.cc> <README.md>")
        sys.exit(1)
    src = Path(sys.argv[1])
    readme = Path(sys.argv[2])
    if not src.exists() or not readme.exists():
        print("Source or README not found")
        sys.exit(1)
    block = generate_block(src)
    ok = replace_in_readme(readme, block)
    if ok:
        print("README updated with function references.")
    else:
        sys.exit(2)


if __name__ == '__main__':
    main()
