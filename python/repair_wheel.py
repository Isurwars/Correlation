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

    if sys.platform == "win32":
        # Find tbb.dll recursively under the build directory
        tbb_dlls = glob.glob("build/**/tbb.dll", recursive=True)
        if not tbb_dlls:
            print("Warning: tbb.dll not found in build directory. Check if it was built statically or not yet compiled.", file=sys.stderr)
            tbb_dir = ""
        else:
            tbb_dir = os.path.dirname(os.path.abspath(tbb_dlls[0]))
            print(f"Found tbb.dll in: {tbb_dir}")

        cmd = ["delvewheel", "repair"]
        if tbb_dir:
            cmd.extend(["--add-path", tbb_dir])
        cmd.extend(["-w", dest_dir, wheel])
        
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    elif sys.platform == "darwin":
        # Find libtbb.dylib recursively under the build directory
        tbb_dylibs = glob.glob("build/**/libtbb*.dylib", recursive=True)
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
        cmd = ["auditwheel", "repair", "-w", dest_dir, wheel]
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
