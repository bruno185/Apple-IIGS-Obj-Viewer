VS 2022 build instructions (GDI viewer)

Options:
- Use CMake to generate a Visual Studio 2022 solution (recommended), or
- Open this folder in Visual Studio (File → Open → Folder) and let VS configure CMake presets.

From "x64 Native Tools Command Prompt for VS 2022":

1) Configure with CMake:
   cmake -S . -B build_vs -G "Visual Studio 17 2022" -A x64

2) Build:
   cmake --build build_vs --config Release --target viewer_win32

3) Run:
   build_vs\bin\Release\viewer_win32.exe q1.obj --h 240

Notes:
- Target `viewer_win32` is the GDI-based native viewer.
- Use Visual Studio's "Open Folder" or the generated solution in `build_vs` to step through and debug.
- Keep SDL-based target (`viewer`) in the repo if you want cross-platform testing.
