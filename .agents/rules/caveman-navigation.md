# Rule: Caveman Navigation Protocol (Token-Efficient Codebase Discovery)

*Activation Mode: Universal / Default Baseline*

## 1. Core Directives

1. **Graphify First (Mandatory Starting Action):**
   - Before opening, reading, or inspecting any source code files (`.cpp`, `.hpp`, `.h`, `CMakeLists.txt`), check `graphify-out/GRAPH_REPORT.md` or `graphify-out/graph.json`.
   - Identify target modules, community clusters, and dependency paths from graph artifacts *before* taking any file-reading action.

2. **Prohibited Baseline Actions:**
   - **No Blind File Reads:** Opening source files directly for speculative discovery or exploration is strictly prohibited.
   - **No Imports Grepping:** Do not grep `#include` directives or class definitions across the workspace to trace dependencies; use `graphify-out/` artifacts instead.

3. **Fallback Protocol (When Graphify Artifacts are Missing/Stale):**
   - If `graphify-out/GRAPH_REPORT.md` does not exist, run `graphify update .` *first* if available, or fall back to high-level directory tree inspection (`ls -R` / `find`) before reading source content.

4. **Targeted Inspection ("Context Diet"):**
   - Treat context tokens as a strict budget. Open *only* the single, minimal set of files identified by the dependency graph required to fulfill the prompt.
   - Do not pull in header files or upstream dependencies unless a compilation or type definition error explicitly demands it.

## 2. Execution Workflow

- **Step 1 (Map):** Consult `graphify-out/GRAPH_REPORT.md` to map the target node or cluster.
- **Step 2 (Filter):** Select the exact file path(s) tied to the relevant node.
- **Step 3 (Execute):** Perform the required analysis or edit on those specific target files only.