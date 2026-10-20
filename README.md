## Building

### Requirements

Install:

- CMake
- Ninja
- GCC via MSYS2 UCRT64

Make sure the compiler is available from PowerShell:

```powershell
g++ --version
cmake --version
ninja --version
```

### Development build

From the project root:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build
```

Run with:

```powershell
.\build\finite_difference\numerical-solver.exe
```

After the initial configuration, rebuilding only requires:

```powershell
cmake --build build
```

### Release build

Configure a separate release build directory:

```powershell
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

Run with:

```powershell
.\build-release\finite_difference\numerical-solver.exe
```

Using separate `build` and `build-release` directories allows the Debug and Release versions to exist at the same time.

### Clean rebuild

Delete the relevant build directory and configure again:

```powershell
Remove-Item -Recurse -Force build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build
```
