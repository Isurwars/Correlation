---
name: cpp-standards
description: Enforce modern C++23 standards, RAII resource management, immutability by default, and memory locality. Use when writing, refactoring, or reviewing modern C++ code.
---

# Modern C++ Standards (C++23)

This skill enforces high-performance, modern C++ (C++23) coding standards across the workspace.

## Core Operational Directives

1. **RAII Management (R.1, R.11):** No naked `new`/`delete` or `malloc`/`free`. Manage all resources via RAII containers (`std::unique_ptr`, `std::shared_ptr`, custom RAII handles).
2. **Const & Immutability (Con.1, Con.2, ES.25):** Default to `const` for variables and `const` for member functions. Use `constexpr` for values computable at compile time.
3. **Type & Conversion Safety (P.4, ES.46, ES.48):** Zero implicit narrowing conversions (`-Werror=float-conversion`, `-Werror=implicit-float-conversion`, `-Werror=shorten-64-to-32`, MSVC `/we4244 /we4267 /we4305`). Wrap `double` literals targeting `real_t` in `static_cast<real_t>(...)`, cast `size_t` to `AtomID` explicitly, use `[](unsigned char chr) { return static_cast<char>(std::tolower(chr)); }` in `std::ranges::transform`, and qualify `std::sqrt`/`std::cos`/`std::sin`/`std::atan2` in C++/CUDA code. Prohibit C-style casts; use `static_cast`, `std::bit_cast`, or explicit constructors.
4. **Parameter Passing (F.16, F.20):** Pass cheap types by value, large objects by `const&`, and sink parameters by value (with `std::move`). Return structs instead of using output parameters.
5. **Class Design (C.20, C.21, C.35, C.46):** Enforce Rule of Zero or Rule of Five. Mark single-argument constructors `explicit`. Base class destructors must be public virtual or protected non-virtual.
6. **Concurrency & OpenMP Safety (CP.2, CP.20, CP.44):** Loop variables in OpenMP constructs must be thread-private. Use `std::scoped_lock` for multi-mutex locking and always name lock guards.
7. **Zero Suppression Policy:** Strictly prohibit `NOLINT`, `NOLINTNEXTLINE`, `NOLINTBEGIN`, and inline suppression comments. Resolve all linter diagnostics structurally.

## Extended Reference

For deep guideline breakdowns, parameter matrices, and anti-pattern examples, consult the internal reference:
- [C++ Core Guidelines Reference](file:///home/isurwars/Projects/Correlation/.agents/skills/cpp-standards/references/core-guidelines.md)
