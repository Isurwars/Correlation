# Skill: cpp-quality-gate
Description: Runs clang-format and clang-tidy over modified files.

## Actions
When code in a C++ file is modified:
1. Run `clang-format -i <filename>` to clean up layout.
2. Run `clang-tidy <filename> -- -Iinclude` (or using `compile_commands.json`) to catch static analysis errors.
3. If errors are found, parse the compiler output, fix the code signature or structure inline, and re-run until it passes clean.