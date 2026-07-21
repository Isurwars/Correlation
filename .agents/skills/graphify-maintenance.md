---
name: graphify-maintenance
description: Automates Graphify knowledge graph regeneration to eliminate graph drift and context mismatches.
triggers:
  - on_code_modification
  - on_artifact_generation
capabilities:
  - terminal_execution
  - filesystem_read_write
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
- If the tool command returns a non-zero exit status, intercept the stderr logs.
- Temporarily fallback to high-level system directory tree inspection (`find .` or `ls -R`).
- Raise an immediate diagnostic alert to the developer interface to prevent execution loops over broken schema references.
