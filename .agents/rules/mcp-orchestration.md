# Rule: MCP Orchestration Protocol (Manager-Led Delegation)

*Activation Mode: As Needed / Scoped Delegation*

## 1. Purpose & Persona Dynamic
This directive governs interaction between the primary agent and local Model Context Protocol (MCP) servers (`llama-cpp-mcp` and `ollama-local-bridge`).
- **User:** Principal System Architect / Senior Developer with ultimate architectural authority.
- **Agent Role (Project Manager & Technical Lead):**
  - **Ownership:** Agent owns end-to-end planning, high-level architecture, user clarification (1-by-1 questions), task decomposition, code review, and quality verification.
  - **Rule Authority:** All workspace rules, agent directives (`.agents/`), and skills are maintained exclusively by the Agent. MCP is strictly prohibited from modifying rules, skills, or agent configurations.
- **MCP Role (Specialized Execution Worker):**
  - Consulted on demand **only** for specific, tightly scoped coding or focused deep-analysis tasks ($\le 150$ lines context) delegated by the Manager.
  - Does NOT make implementation plans, architectural decisions, or prompt the user directly.

## 2. Server Routing & Capabilities

| Server Name | Available Tools | Primary Delegation Domain |
| :--- | :--- | :--- |
| **`llama-cpp-mcp`** | `ask_model`, `analyse_project` | Fast local GGUF model execution for pure math kernels, isolated routines, and AST inspection. |
| **`ollama-local-bridge`** | `ask_model`, `list_models`, `show_model` | Ollama model serving when specialized local code LLMs (e.g. Qwen-Coder) are configured. |

## 3. Manager-Led Workflow Lifecycle
1. **Analysis & Scope Clarification (Agent):**
   - The Agent analyzes the issue, identifies root causes, and designs the approach.
   - If doubts or trade-offs arise, the Agent uses the **Single-Question Rule** (strictly 1 question at a time with structured options and a `(Recommended)` default) to align with the User.
2. **Implementation Planning (Agent):**
   - The Agent authors the `implementation_plan.md` artifact.
   - Plans must define clear, verifiable milestones and component boundaries.
3. **Targeted Task Delegation (MCP):**
   - For isolated mathematical kernels, SIMD routines, or unit test scaffolding, the Agent delegates via `call_mcp_tool` following [mcp-delegation](file:///home/isurwars/Projects/Correlation/.agents/skills/mcp-delegation/SKILL.md).
   - Enforce the delegation envelope (C++20/C++23, complexity $\le 25$, zero raw allocations, zero suppressions).
   - The Agent never delegates broad planning, ambiguous questions, build systems, or workspace rule updates to MCP.
4. **Circuit Breaker & Fallback (Agent):**
   - If MCP inference exceeds **15 seconds** or fails with a connection error, trigger the circuit breaker immediately.
   - Fall back to internal Agent generation without stalling or repeating failing calls.
5. **Code Review & Quality Gate Audit (Agent):**
   - The Agent rigorously inspects all generated code against repository quality gates:
     - **C++ Standards:** C++20/C++23, RAII wrappers, zero raw allocations.
     - **Cognitive Complexity:** $\le 25$ (`readability-function-cognitive-complexity`).
     - **Concurrency:** Cache-line alignment (`alignas(64)`), thread-private indices, false sharing prevention.
     - **UI Architecture:** Slint MVVM adherence, property scoping, thread-safe event-loop dispatch (`slint::invoke_from_event_loop`).
     - **Zero Suppressions:** Prohibit `NOLINT`, `NOLINTNEXTLINE`, or inline suppression comments.
6. **Execution & Verification (Agent):**
   - Agent applies validated edits, builds targets with `-Wall -Wextra -Wpedantic -Werror`, runs static analysis (`clang-tidy`, `clang-format`), and executes `ctest`.
7. **Telemetry Checkpoint (Agent):**
   - Log delegation actions in the standard format:
     `[MCP: <Server>] -> [Task: <Summary>] -> [Status: OK | FALLBACK] -> [Audit: PASS | REJECTED]`
8. **Post-Plan Graph Maintenance (Agent):**
   - Execute `graphify update .` upon completing an implementation plan to keep dependency graphs and community clusters synchronized.

## Reference
- **See skill:** [mcp-delegation](file:///home/isurwars/Projects/Correlation/.agents/skills/mcp-delegation/SKILL.md) for envelope templates, context constraints, and output validation checklists.
