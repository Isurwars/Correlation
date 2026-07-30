# -*- coding: utf-8 -*-
# @file repair_wheel.py
# @brief Helper script to dynamically resolve and vendor oneTBB dependencies for wheels.
# @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# @par License
# SPDX-License-Identifier: AGPL-3.0-only

import os
import sys
import glob
import subprocess

def main():
    if len(sys.argv) < 3:
        print("Usage: repair_wheel.py <wheel_path> <dest_dir>", file=sys.stderr)
        sys.exit(1)

    wheel = sys.argv[1]
    dest_dir = sys.argv[2]

    print(f"Repairing wheel: {wheel}")
    print(f"Destination: {dest_dir}")

    wheel_name = os.path.basename(wheel).lower()

    if sys.platform == "win32":
        # Find tbb DLLs recursively under the build directory matching target arch
        tbb_dlls = glob.glob("build/**/*tbb*.dll", recursive=True)
        matching = [f for f in tbb_dlls if ("amd64" in f.lower() or "x64" in f.lower() or "win64" in f.lower())] if "amd64" in wheel_name else tbb_dlls
        if matching:
            tbb_dlls = matching

        if not tbb_dlls:
            print("Warning: tbb DLLs not found in build directory. Check if built statically or not yet compiled.", file=sys.stderr)
            tbb_dirs = []
        else:
            tbb_dirs = sorted(list({os.path.dirname(os.path.abspath(f)) for f in tbb_dlls}))
            print(f"Found tbb DLL directories: {tbb_dirs}")

        cmd = ["delvewheel", "repair"]
        for d in tbb_dirs:
            cmd.extend(["--add-path", d])
        cmd.extend(["-w", dest_dir, wheel])
        
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    elif sys.platform == "darwin":
        # Find libtbb.dylib recursively under the build directory matching target arch
        tbb_dylibs = glob.glob("build/**/libtbb*.dylib", recursive=True)
        if "arm64" in wheel_name:
            matching = [f for f in tbb_dylibs if "arm64" in f]
        elif "x86_64" in wheel_name:
            matching = [f for f in tbb_dylibs if "x86_64" in f]
        else:
            matching = []

        if matching:
            tbb_dylibs = matching

        if not tbb_dylibs:
            print("Warning: libtbb.dylib not found in build directory.", file=sys.stderr)
            tbb_dir = ""
        else:
            tbb_dir = os.path.dirname(os.path.abspath(tbb_dylibs[0]))
            print(f"Found libtbb.dylib in: {tbb_dir}")

        env = os.environ.copy()
        if tbb_dir:
            current_dyld = env.get("DYLD_LIBRARY_PATH", "")
            env["DYLD_LIBRARY_PATH"] = f"{tbb_dir}:{current_dyld}" if current_dyld else tbb_dir
            print(f"Setting DYLD_LIBRARY_PATH={env['DYLD_LIBRARY_PATH']}")

        cmd = ["delocate-wheel", "-w", dest_dir, wheel]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, env=env, check=True)

    else:
        # Linux (auditwheel)
        tbb_sos = glob.glob("build/**/libtbb*.so*", recursive=True)
        if "aarch64" in wheel_name:
            matching = [f for f in tbb_sos if "aarch64" in f]
        elif "x86_64" in wheel_name:
            matching = [f for f in tbb_sos if "x86_64" in f]
        else:
            matching = []

        if matching:
            tbb_sos = matching

        if not tbb_sos:
            print("Warning: libtbb.so not found in build directory.", file=sys.stderr)
            tbb_dir = ""
        else:
            tbb_dirs = sorted(list({os.path.dirname(os.path.abspath(f)) for f in tbb_sos}))
            tbb_dir = ":".join(tbb_dirs)
            print(f"Found libtbb shared libraries in: {tbb_dir}")

        env = os.environ.copy()
        if tbb_dir:
            current_ld = env.get("LD_LIBRARY_PATH", "")
            env["LD_LIBRARY_PATH"] = f"{tbb_dir}:{current_ld}" if current_ld else tbb_dir
            print(f"Setting LD_LIBRARY_PATH={env['LD_LIBRARY_PATH']}")

        cmd = ["auditwheel", "repair", "-w", dest_dir, wheel]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, env=env, check=True)

if __name__ == "__main__":
    main()
