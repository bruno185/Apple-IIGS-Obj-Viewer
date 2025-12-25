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

Build (Windows, using CMake + Visual Studio or MSYS2)

There are convenience scripts provided in this directory:
- `build_msvc.ps1` — try to generate Visual Studio project and build (uses CMake generators for Visual Studio 2019/2022)
- `build_mingw.ps1` — generates MinGW Makefiles and builds
- `smoke_test.ps1` — runs a built executable with `test_small.obj` and checks the log for startup messages

Manual steps (CMake):
1. Install CMake and a Windows toolchain (MSVC or MinGW)
2. Choose generator and build directory, for example:
   ```powershell
   mkdir build_mingw
   cmake -S . -B build_mingw -G "MinGW Makefiles"
   cmake --build build_mingw --config Release
   ```
3. Run the built exe (example):
   ```powershell
   .\build_mingw\bin\viewer_win32.exe ..\..\test_small.obj
   ```

If you prefer Visual Studio: use `build_msvc.ps1` or run `cmake -S . -B build_vs -G "Visual Studio 17 2022" -A x64` and open the generated solution.

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

