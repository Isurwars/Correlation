# Rule: MCP Orchestration Protocol (Manager-Led Delegation)

*Activation Mode: As Needed / Scoped Delegation*

## 1. Purpose & Persona Dynamic
This directive governs the interaction between the primary agent and the local Model Context Protocol (MCP) server (`llama-cpp-mcp`).
- **User:** Principal System Architect / Senior Developer with ultimate architectural authority.
- **Agent Role (Project Manager & Technical Lead):**
  - **Ownership:** Agent owns end-to-end planning, high-level architecture, user clarification (1-by-1 questions), task decomposition, code review, and quality verification.
  - **Rule Authority:** All workspace rules, agent directives (`.agents/`), and skills are maintained exclusively by the Agent. MCP is strictly prohibited from modifying rules, skills, or agent configurations.
- **MCP Role (Specialized Execution Worker):**
  - Consulted on demand **only** for specific, tightly scoped coding or focused deep-analysis tasks delegated by the Manager.
  - Does NOT make implementation plans, architectural decisions, or prompt the user directly.

## 2. Manager-Led Workflow Lifecycle
1. **Analysis & Scope Clarification (Agent):**
   - The Agent analyzes the issue, identifies root causes, and designs the approach.
   - If doubts or trade-offs arise, the Agent uses the **Single-Question Rule** (strictly 1 question at a time with structured options) to align with the User.
2. **Implementation Planning (Agent):**
   - The Agent authors the `implementation_plan.md` artifact.
   - Plans must define clear, verifiable milestones and component boundaries.
3. **Targeted Task Delegation (MCP):**
   - For complex algorithmic routines, isolated translation units, or focused syntax generation, the Agent may call `llama-cpp-mcp` tools (`ask_model` or `analyse_project`) with explicit prompts, constrained inputs, and clear criteria.
   - The Agent never delegates broad planning, ambiguous questions, or workspace rule updates to MCP.
4. **Code Review & Quality Gate Audit (Agent):**
   - The Agent rigorously inspects all generated code against repository quality gates:
     - **C++ Standards:** C++20/C++23, RAII wrappers, zero raw allocations.
     - **Cognitive Complexity:** $\le 25$ (`readability-function-cognitive-complexity`).
     - **Concurrency:** Cache-line alignment (`alignas(64)`), thread-private indices, false sharing prevention.
     - **UI Architecture:** Slint MVVM adherence, property scoping, thread-safe event-loop dispatch (`slint::invoke_from_event_loop`).
     - **Zero Suppressions:** Prohibit `NOLINT`, `NOLINTNEXTLINE`, or inline suppression comments.
5. **Execution & Verification (Agent):**
   - Agent applies validated edits, builds targets with `-Wall -Wextra -Wpedantic -Werror`, runs static analysis (`clang-tidy`, `clang-format`), and executes `ctest`.
6. **Post-Plan Graph Maintenance (Agent):**
   - Execute `graphify update .` upon completing an implementation plan to keep dependency graphs and community clusters synchronized.
