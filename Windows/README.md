Windows viewer for GS3Dp algorithm

Overview
--------
This is a compact C implementation (SDL2) that reproduces the GS3Dp pairwise ordering algorithm and renders polygons in the computed painter order on Windows.

Features
- OBJ loader (v and f simple format)
- Fixed32 16.16 helpers to keep behaviour close to the original
- Pairwise comparator approximating GS3Dp logic (depth, bbox, plane votes, triangle votes, simple local sampling)
- SDL2-based renderer that displays filled polygons in order
- Diagnostic CSV (PAIRDIAG.CSV) compatible lines

Build (Windows, recommended: MSVC via included script)

Quick (recommended):
- From the repository root, run the convenience script that generates and builds a Release MSVC solution:

```powershell
powershell -ExecutionPolicy Bypass -File "Windows\build_msvc.ps1"
```

This script configures CMake (tries Visual Studio 2019/2022) and builds the `viewer_win32` Release target.

Manual (CMake + Visual Studio):
1. Ensure CMake and a Visual Studio toolchain are installed (e.g. Visual Studio 2022).
2. Configure and build:

```powershell
cmake -S Windows -B Windows/build_vs -G "Visual Studio 17 2022" -A x64
cmake --build Windows/build_vs --config Release --target viewer_win32
```

3. Run the executable:

```powershell
Windows\build_vs\bin\Release\viewer_win32.exe <path-to-obj> [--angle_h <deg>] [--angle_v <deg>] [--distance <dist>]
# example:
Windows\build_vs\bin\Release\viewer_win32.exe q1.obj
```

MinGW: there is no maintained `build_mingw.ps1` in this repository; if you need a MinGW build, generate MinGW Makefiles with CMake and build them manually (see `cmake -G "MinGW Makefiles"`).

Smoke test
----------
- Use `Windows\smoke_test.ps1` to run a simple startup test; it checks that logs are written to `%TEMP%\viewer_win32.log` and that the viewer starts successfully.

Notes
-----
- Logs and runtime diagnostics are written to `%TEMP%\viewer_win32.log`.
- The Windows target `viewer_win32` is a GDI-based native viewer and is the primary Windows target.
- If you prefer to open the generated solution, use VS "Open Project/Solution" on `Windows/build_vs\GS3DpViewer.sln` (after configuring with CMake).

Smoke test
----------
Run `.uild_mingw.ps1` (or `.uild_msvc.ps1`) and then `.\nsmoke_test.ps1` to verify the executable starts and writes diagnostic logs to `%TEMP%\viewer_win32.log`.

Notes
-----
- The Win32 frontend does not use ORCA — it is a native Windows build and is intended to be built on Windows.
- The painter / face ordering code is preserved in `src/painter.c` to maintain algorithm consistency with the original project.


Running
-------
viewer.exe <path-to-obj> [--h <height>] [--angle_h <deg>] [--angle_v <deg>] [--distance <dist>]

Notes
-----
This is an initial faithful port focusing on correctness and diagnostics. If you need bit-identical behaviour with GS3Dp, we can further tune fixed-point math and epsilons.

Files
-----
- CMakeLists.txt
- src/main.c
- src/fixed32.h
- src/obj.c
- src/painter.c
- src/renderer.c

