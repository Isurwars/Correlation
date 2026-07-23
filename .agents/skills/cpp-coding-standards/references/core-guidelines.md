# C++ Core Guidelines Reference Manual

Comprehensive coding standards for modern C++ (C++17/20/23) derived from the [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines). Enforces type safety, resource safety, immutability, and clarity.

## Cross-Cutting Principles

1. **RAII everywhere** (P.8, R.1, E.6, CP.20): Bind resource lifetime to object lifetime
2. **Immutability by default** (P.10, Con.1-5, ES.25): Start with `const`/`constexpr`; mutability is the exception
3. **Type safety** (P.4, I.4, ES.46-49, Enum.3): Use the type system to prevent errors at compile time
4. **Express intent** (P.3, F.1, NL.1-2, T.10): Names, types, and concepts should communicate purpose
5. **Minimize complexity** (F.2-3, ES.5, Per.4-5): Simple code is correct code
6. **Value semantics over pointer semantics** (C.10, R.3-5, F.20, CP.31): Prefer returning by value and scoped objects

---

## Philosophy & Interfaces (P.*, I.*)

### Key Rules

| Rule     | Summary                                                |
| -------- | ------------------------------------------------------ |
| **P.1**  | Express ideas directly in code                         |
| **P.3**  | Express intent                                         |
| **P.4**  | Ideally, a program should be statically type safe      |
| **P.5**  | Prefer compile-time checking to run-time checking      |
| **P.8**  | Don't leak any resources                               |
| **P.10** | Prefer immutable data to mutable data                  |
| **I.1**  | Make interfaces explicit                               |
| **I.2**  | Avoid non-const global variables                       |
| **I.4**  | Make interfaces precisely and strongly typed           |
| **I.11** | Never transfer ownership by a raw pointer or reference |
| **I.23** | Keep the number of function arguments low              |

```cpp
// P.10 + I.4: Immutable, strongly typed interface
struct Temperature {
    double kelvin;
};

Temperature boil(const Temperature& water);
```

---

## Functions (F.*)

### Key Rules

| Rule     | Summary                                                                        |
| -------- | ------------------------------------------------------------------------------ |
| **F.1**  | Package meaningful operations as carefully named functions                     |
| **F.2**  | A function should perform a single logical operation                           |
| **F.3**  | Keep functions short and simple                                                |
| **F.4**  | If a function might be evaluated at compile time, declare it `constexpr`       |
| **F.6**  | If your function must not throw, declare it `noexcept`                         |
| **F.8**  | Prefer pure functions                                                          |
| **F.16** | For "in" parameters, pass cheaply-copied types by value and others by `const&` |
| **F.20** | For "out" values, prefer return values to output parameters                    |
| **F.21** | To return multiple "out" values, prefer returning a struct                     |
| **F.43** | Never return a pointer or reference to a local object                          |

---

## Classes & Class Hierarchies (C.*)

### Key Rules

| Rule      | Summary                                                                             |
| --------- | ----------------------------------------------------------------------------------- |
| **C.2**   | Use `class` if invariant exists; `struct` if data members vary independently        |
| **C.9**   | Minimize exposure of members                                                        |
| **C.20**  | If you can avoid defining default operations, do (Rule of Zero)                     |
| **C.21**  | If you define or `=delete` any copy/move/destructor, handle them all (Rule of Five) |
| **C.35**  | Base class destructor: public virtual or protected non-virtual                      |
| **C.41**  | A constructor should create a fully initialized object                              |
| **C.46**  | Declare single-argument constructors `explicit`                                     |
| **C.67**  | A polymorphic class should suppress public copy/move                                |
| **C.128** | Virtual functions: specify exactly one of `virtual`, `override`, or `final`         |

---

## Resource Management (R.*)

### Key Rules

| Rule     | Summary                                                        |
| -------- | -------------------------------------------------------------- |
| **R.1**  | Manage resources automatically using RAII                      |
| **R.3**  | A raw pointer (`T*`) is non-owning                             |
| **R.5**  | Prefer scoped objects; don't heap-allocate unnecessarily       |
| **R.10** | Avoid `malloc()`/`free()`                                      |
| **R.11** | Avoid calling `new` and `delete` explicitly                    |
| **R.20** | Use `unique_ptr` or `shared_ptr` to represent ownership        |
| **R.21** | Prefer `unique_ptr` over `shared_ptr` unless sharing ownership |
| **R.22** | Use `make_shared()` to make `shared_ptr`s                      |

---

## Expressions & Statements (ES.*)

### Key Rules

| Rule      | Summary                                                                |
| --------- | ---------------------------------------------------------------------- |
| **ES.5**  | Keep scopes small                                                      |
| **ES.20** | Always initialize an object                                            |
| **ES.23** | Prefer `{}` initializer syntax                                         |
| **ES.25** | Declare objects `const` or `constexpr` unless modification is intended |
| **ES.28** | Use lambdas for complex initialization of `const` variables            |
| **ES.45** | Avoid magic constants; use symbolic constants                          |
| **ES.46** | Avoid narrowing/lossy arithmetic conversions                           |
| **ES.47** | Use `nullptr` rather than `0` or `NULL`                                |
| **ES.48** | Avoid casts                                                            |
| **ES.50** | Don't cast away `const`                                                |

---

## Error Handling (E.*)

### Key Rules

| Rule     | Summary                                                                      |
| -------- | ---------------------------------------------------------------------------- |
| **E.1**  | Develop an error-handling strategy early in a design                         |
| **E.2**  | Throw an exception to signal that a function can't perform its assigned task |
| **E.6**  | Use RAII to prevent leaks                                                    |
| **E.12** | Use `noexcept` when throwing is impossible or unacceptable                   |
| **E.14** | Use purpose-designed user-defined types as exceptions                        |
| **E.15** | Throw by value, catch by reference                                           |
| **E.16** | Destructors, deallocation, and swap must never fail                          |
| **E.17** | Don't try to catch every exception in every function                         |

---

## Concurrency & Parallelism (CP.*)

### Key Rules

| Rule       | Summary                                                       |
| ---------- | ------------------------------------------------------------- |
| **CP.2**   | Avoid data races                                              |
| **CP.3**   | Minimize explicit sharing of writable data                    |
| **CP.4**   | Think in terms of tasks, rather than threads                  |
| **CP.8**   | Don't use `volatile` for synchronization                      |
| **CP.20**  | Use RAII, never plain `lock()`/`unlock()`                     |
| **CP.21**  | Use `std::scoped_lock` to acquire multiple mutexes            |
| **CP.22**  | Never call unknown code while holding a lock                  |
| **CP.42**  | Don't wait without a condition                                |
| **CP.44**  | Remember to name your `lock_guard`s and `unique_lock`s        |
| **CP.100** | Don't use lock-free programming unless you absolutely have to |

---

## Templates & Generic Programming (T.*)

### Key Rules

| Rule      | Summary                                                     |
| --------- | ----------------------------------------------------------- |
| **T.1**   | Use templates to raise the level of abstraction             |
| **T.2**   | Use templates to express algorithms for many argument types |
| **T.10**  | Specify concepts for all template arguments                 |
| **T.11**  | Use standard concepts whenever possible                     |
| **T.13**  | Prefer shorthand notation for simple concepts               |
| **T.43**  | Prefer `using` over `typedef`                               |
| **T.120** | Use template metaprogramming only when you really need to   |
| **T.144** | Don't specialize function templates (overload instead)      |
