---
name: graphify-maintenance
description: Automates Graphify knowledge graph regeneration to eliminate graph drift and context mismatches.
---

# Graphify Automation & Sync Protocol

This skill maintains the structural dependency graph (`graphify-out/`) for the Correlation codebase, ensuring agent context stays synchronized with source code changes.

## 1. Trigger Conditions

Run graph maintenance whenever any of the following files are added, modified, or deleted:
- C++ source files (`*.cpp`, `*.hpp`, `*.h`)
- Slint UI components (`*.slint`)
- Build scripts (`CMakeLists.txt`, `*.cmake`)

## 2. Execution Protocol

### Step 1: Run Graphify Update
Execute the codebase graph update command from the project root:
```bash
graphify update .
```

### Step 2: Artifact Validation
Confirm the following output artifacts exist and are non-empty:
1. `graphify-out/graph.json`: Machine-readable dependency graph (must contain >4000 nodes).
2. `graphify-out/GRAPH_REPORT.md`: Human-readable community cluster report (must contain navigation hubs).

### Step 3: Stale Node Verification
If files were deleted or renamed, verify that obsolete nodes no longer appear in `graphify-out/GRAPH_REPORT.md`.

## 3. Failure Handling & Fallbacks

If `graphify` is unavailable or returns an error:
1. Log the failure reason.
2. Fall back to manual target mapping using `CMakeLists.txt` and file tree inspection (`find src/ include/ ui/ -type f`).
3. Alert the user if critical dependency maps are out of sync.
