# 3D OBJ Viewer (Apple IIGS)

A 3D viewer for the Apple IIGS written in ORCA/C using Fixed32 (16.16) arithmetic.

## Description
This project reads simplified OBJ files (vertices `v` and faces `f`), performs 3D transformations, projects to 2D with perspective, and draws polygons using the Apple IIGS QuickDraw API.

## Main features
- Read OBJ files (vertices and faces)
- Optimized 3D transforms using Fixed32 (16.16)
- Fast trig tables for speed
- Painter's algorithm with Newell/Sancha overlap tests
- Non-destructive auto-scale on import (optional) with revert (press `r`)
- Interactive options (angles, distance, color palette)

## Usage
1. Compile with your ORCA/C toolchain or a compatible tool (examples in the repo). Example:

   iix compile GS3Dp.cc

2. Run the generated executable. The program will ask for an OBJ filename and optionally apply auto-scale.

3. Keyboard controls:
- Space: show parameters and redraw (shows whether auto-scale is ON and its factor)
- A / Z: decrease/increase distance (10% steps)
- +/-: apply auto-fit if none present then increase/decrease distance (10% steps)
- Arrow keys: adjust angles
- W / X: rotate screen
- C: change palette
- F: toggle frame-only polygons (default: OFF — polygons are filled by default)
- R: revert auto-scale (if applied)
- K: edit angles/distance interactively without reloading the model (ENTER may trigger auto-fit)

## Implementation notes
- Critical computations are optimized to reduce floating conversions and avoid overflow (heavy use of `Fixed32` and `Fixed64`).
- An orientation fix was added for OBJ Z-up exports (swap Y/Z at import) and can be reverted manually.
- Auto‑fit uses a precomputed **bounding sphere** (centroid + radius) for O(1) distance estimates via `computeDistanceFromBoundingSphere()`; `computeDistanceToFit()` (per‑vertex scan) is retained but **deprecated** as a fallback.

## Credits
- Main author: Bruno
- Tribute: *A tribute to Robert DONY* — Author of "Calcul des parties cachées" (Masson, 1986)

---

If you'd like a more detailed build guide, tests, or a license added, say so and I'll add it.