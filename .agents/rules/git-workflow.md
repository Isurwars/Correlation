# Rule: Git Workflow & Conventional Commits

*Activation Mode: Universal / Default Baseline*

## 1. Branching Model (Trunk-Based Development)

- **Primary Branch:** `main` is the single source of truth. All production code lives here.
- **Direct Commits:** Small, self-contained changes (bug fixes, documentation, minor refactors) commit directly to `main` after local test verification.
- **Feature Branches:** Use short-lived branches (`feat/<name>`, `fix/<name>`) only for multi-day experimental work or breaking changes. Merge back to `main` within 1–3 days.
- **Release Tagging:** Tag releases on `main` using semantic versioning: `git tag -a v1.2.3 -m "Release v1.2.3"`.

## 2. Commit Message Format (Conventional Commits)

All commit messages must follow the [Conventional Commits](https://www.conventionalcommits.org/) specification:

```
<type>(<scope>): <description>

[optional body: explain WHY, not WHAT]

[optional footer: BREAKING CHANGE: <description>]
```

### Permitted Types

| Type | Usage |
| :--- | :--- |
| `feat` | New feature or capability |
| `fix` | Bug fix |
| `refactor` | Code restructuring without behavior change |
| `perf` | Performance improvement |
| `test` | Adding or updating tests |
| `docs` | Documentation only |
| `build` | Build system or dependency changes (CMake, FetchContent) |
| `ci` | CI/CD pipeline changes |
| `chore` | Maintenance tasks (formatting, linting, agent hooks) |

### Permitted Scopes

- **Calculators / Physics:** `rdf`, `pdf`, `sq`, `pad`, `rings`, `msd`, `vacf`, `vdos`, `lef`, `calculators`, `analysis`, `math`, `physics`, `mlip`
- **I/O & Formats:** `readers`, `writers`, `io`
- **UI / Presentation:** `ui`, `slint`, `plotters`
- **Interfaces & Infrastructure:** `bindings`, `cli`, `core`, `utils`, `cmake`, `tests`, `docs`

### Breaking Changes
Append `!` after the type: `feat(core)!: rename DistributionFunctions to CorrelationCalculator`

## 3. Prohibited Actions
- **Never** force-push to `main`.
- **Never** commit generated files (`build/`, `build-*/`, `*.o`, `compile_commands.json`, `graphify-out/`).
- **Never** commit secrets, API keys, or credentials.
