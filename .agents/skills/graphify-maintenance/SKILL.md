---
name: graphify-maintenance
description: Automates Graphify knowledge graph regeneration to eliminate graph drift and context mismatches.
---

# Graphify Automation & Sync Protocol

## 1. Core Update Trigger
- **Condition**: Every time a source file (`.cpp`, `.hpp`, `.h`, `CMakeLists.txt`) is successfully written, modified, or patched by an agent.
- **Action**: Immediately execute the codebase graph update command to synchronize dependencies before any subsequent task begins.

## 2. Execution Blueprint

### Step 1: Run Update Command
Execute the graphify update command relative to the project root directory:
```bash
graphify update .
```

### Step 2: Validate Outputs
Verify that the underlying file system artifacts are fresh:
- Check that `graphify-out/graph.json` has a modified timestamp matching the current system runtime execution window.
- Verify that `graphify-out/GRAPH_REPORT.md` is populated with correct, updated structural nodes.

## 3. Handling Update Failures
- If `graphify` command is missing or returns a non-zero exit status, intercept the logs.
- Fall back to directory tree inspection (`find .` or `ls -R`).
- Raise a diagnostic alert if schema references are broken.
