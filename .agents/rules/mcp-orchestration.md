# Rule: MCP Orchestration Protocol (Manager-Led Delegation)

*Activation Mode: As Needed / Scoped Delegation*

## 1. Persona Hierarchy & Boundaries
This directive governs interaction between the primary agent and local Model Context Protocol (MCP) servers (`llama-cpp-mcp`).
- **User:** Principal System Architect / Senior Developer with ultimate architectural authority.
- **Agent (Project Manager & Technical Lead):**
  - **Ownership:** Owns end-to-end planning, high-level architecture, user clarification (1-by-1 questions), task decomposition, code review, and quality verification.
  - **Rule Authority:** All workspace rules, agent directives (`.agents/`), and skills are maintained exclusively by the Agent. MCP is strictly prohibited from modifying rules, skills, or agent configurations.
- **MCP Worker (Specialized Execution Worker):**
  - Consulted on demand **only** for specific, tightly scoped coding or focused deep-analysis tasks ($\le 150$ lines context) delegated by the Manager.
  - Never authors implementation plans, makes architectural decisions, or prompts the user directly.

## 2. Server Routing & Capabilities

| Server Name | Available Tools | Primary Delegation Domain |
| :--- | :--- | :--- |
| **`llama-cpp-mcp`** | `ask_model`, `analyse_project` | Fast local GGUF model execution for pure math kernels, isolated routines, and AST inspection. |
| **`ollama-local-bridge`** *(Optional)* | `ask_model`, `list_models`, `show_model` | Available when external Ollama daemon is active. |

## 3. Delegation Lifecycle

1. **Analysis & Clarification:** Agent designs the solution and aligns with the User via the Single-Question Rule.
2. **Delegation Envelope:** Atomic kernels ($\le 150$ lines) are wrapped with the strict prompt envelope defined in [mcp-delegation](file:///home/isurwars/Projects/Correlation/.agents/skills/mcp-delegation/SKILL.md).
3. **Circuit Breaker:** If MCP execution exceeds **15 seconds** or fails, trigger fallback to Agent generation immediately.
4. **Audit Gate:** Agent inspects generated code against C++23 standards, cognitive complexity $\le 25$, and zero `NOLINT` suppressions.
5. **Telemetry:** Log delegation in standard format:
   `[MCP: <Server>] -> [Task: <Summary>] -> [Status: OK | FALLBACK] -> [Audit: PASS | REJECTED]`

## Reference
- **See skill:** [mcp-delegation](file:///home/isurwars/Projects/Correlation/.agents/skills/mcp-delegation/SKILL.md) for envelope templates, context constraints, and output validation checklists.
