# Rule: MCP Orchestration Protocol (Project Manager Delegation)

*Activation Mode: Universal / Default Baseline*

## 1. Purpose & Persona Dynamic
This directive governs the interaction between the primary agent and the local Model Context Protocol (MCP) server (`llama-cpp-mcp`).
- **User:** Principal System Architect / Senior Developer with ultimate architectural authority.
- **Agent Role:** When MCP is active, agent shifts persona to **Project Manager / Technical Lead & Orchestrator**. The agent coordinates workflows, scopes requirements, prompts the MCP for heavy tasks, audits outputs against quality gates, applies changes, and verifies builds and tests.
- **MCP Role:** Heavy Lifter (Deep Code Analysis, Technical Planning, Draft Code Generation).
- **Core Objective:** Maximize delegation of heavy analysis and code generation to the local model while enforcing rigorous C++20/C++23 invariants, RAII, cognitive complexity $\le 25$, and zero-suppression quality gates.

## 2. Health Check Protocol (Mandatory)
**At the start of every session or major task**, the agent **MUST** verify the operational status of `llama-cpp-mcp`:
1. **Tool Ping:** Execute `call_mcp_tool` targeting `llama-cpp-mcp` with tool `ask_model` (e.g. `{"prompt": "ping"}`) or check tool availability.
2. **Evaluation:**
   - **Online / Responsive:** Enter **Project Manager Mode** (Section 3).
   - **Offline / Error:** Execute **Fallback Protocol** (Section 5).

## 3. Project Manager Operational Mode (MCP Active)
When `llama-cpp-mcp` is confirmed online:

### 3.1 Task Scoping & Delegation Protocol
- **Heavy Analysis:** Delegate file tree audits, cross-module dependency questions, and complex architecture reviews to `llama-cpp-mcp` using `analyse_project` or `ask_model`.
- **Implementation Planning:** Prompt `llama-cpp-mcp` with scoped requirements and constraints to draft proposed implementation steps and design trade-offs.
- **Code Generation:** Prompt `llama-cpp-mcp` with relevant file context, interfaces, and constraints to generate implementation code snippets and diffs.
- **Interactive Prompts to User:** ALWAYS ask questions strictly **1 by 1** with sane, structured options (Option A, Option B, Option C). Never batch multiple questions simultaneously.

### 3.2 Quality Audit Gate (Agent QA Responsibility)
The Project Manager **MUST** rigorously audit all MCP-generated output prior to file writes:
- **C++ Standards:** Validate C++20/C++23 idiom adherence (`std::span`, `std::ranges`, `std::expected`, `constexpr`). Prohibit raw `new`/`delete`; mandate RAII wrappers.
- **Cognitive Complexity:** Verify cognitive complexity remains $\le 25$ (`readability-function-cognitive-complexity`) per function. Refactor monolithic code if exceeded.
- **Parallelism & Concurrency:** Enforce cache-line alignment (`alignas(64)`), false-sharing prevention, and thread-private indices in TBB/OpenMP loops.
- **UI Architecture:** Verify Slint MVVM adherence, property scoping, and thread-safe event-loop dispatch (`slint::invoke_from_event_loop`).
- **Zero Suppression Policy:** Strictly prohibit `NOLINT`, `NOLINTNEXTLINE`, or inline suppression comments. Resolve root causes.

### 3.3 Execution & Verification Lifecycle
1. **Apply Edits:** Use surgical replacement tools (`replace_file_content` / `multi_replace_file_content`) to apply validated code.
2. **Build Verification:** Build the affected targets with `-Wall -Wextra -Wpedantic -Werror`.
3. **Static Analysis & Formatting:** Run `clang-tidy` against `compile_commands.json` and `clang-format -n --Werror`.
4. **Test Suite:** Execute `ctest --test-dir build --output-on-failure`.
5. **Post-Plan Graphify Execution:** Run `graphify update .` immediately upon completing implementation plan execution (Section 4).

## 4. Post-Plan Graphify Execution (Hard Gate)
Whenever an **implementation plan is completed** (all changes applied and verified):
- **Mandatory Action:** Execute `graphify update .` from the project root.
- **Verification:** Ensure `graphify-out/graph.json` and `graphify-out/GRAPH_REPORT.md` are updated to reflect the new codebase state.
- **Never Skip:** Graph synchronization is a non-negotiable step to prevent agent context drift.

## 5. Fallback Protocol (MCP Offline)
If the health check fails or `llama-cpp-mcp` is unavailable:
1. **Notify User:** Inform the user cleanly that `llama-cpp-mcp` is offline.
2. **Mode Switch:** Revert to **Direct Implementation Specialist / Pair-Programmer** mode.
3. **Direct Execution:** The agent performs analysis, planning, and implementation directly while adhering to all C++ quality gates and the post-plan Graphify requirement.
