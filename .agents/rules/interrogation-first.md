# Rule: Interrogation-First Protocol (Mandatory Pre-Implementation Gate)

*Activation Mode: Universal / Default Baseline*

## 1. Non-Negotiable Constraint

Before writing, modifying, or generating ANY source code, build files, or production artifacts, the agent MUST interrogate the user to extract precise technical requirements.

**Zero Premature Implementation:** Strictly forbidden from modifying source code, creating production files, or executing build actions until requirements are validated through interrogation.

## 2. Single-Question Constraint (Strict 1-by-1 Protocol)

Ask strictly **ONE question at a time**. Never batch multiple questions, combine multiple topics, or present menus of inquiries simultaneously.

- **Sane, Structured Options:** Present 2–3 concrete, well-grounded options (Option A, Option B, Option C) with explicit technical trade-offs (latency, cache locality, memory overhead, API ergonomics). Avoid nonsensical, ungrounded, or frivolous options.
- **Interactive Tooling:** ALWAYS use the `ask_question` tool when soliciting multiple-choice responses from the user, presenting exactly one question per modal.
- **Wait State:** Emit the single question and immediately wait for the user's response before formulating any subsequent inquiries or actions.
- **Proactive Coverage:** Cover non-functional requirements proactively (performance bounds, scale limits, thread safety, UI state synchronization) one topic per turn.

## 3. Proactive Architectural Auditing

Audit every request for blind spots before proposing solutions:
- **Edge cases:** Out-of-bounds inputs, empty vectors/histograms, numerical underflow/overflow.
- **Concurrency hazards:** Race conditions, mutex lock scope, atomic memory ordering, thread-local accumulation.
- **Memory lifespans:** RAII wrappers, zero-copy buffers, ownership transfer.
- **UI state synchronization:** Slint event loop dispatch, VectorModel invalidation, renderer thread safety.

## 4. Understanding Lock (Hard Gate)

Before proposing any implementation plan or code modification, present an **Understanding Lock** summary:

- **What is being built**: Concise definition.
- **Why it exists**: Core motivation.
- **Target users / callers**: UI, C++ core, Python bindings, CLI.
- **Key constraints**: Threading, memory, performance, backwards compatibility.
- **Explicit non-goals**: Scope boundaries.

### Hard Gate Confirmation Prompt:
> *"Does this accurately reflect your intent? Please confirm or correct anything before we finalize the design."*

**Do NOT proceed to design or implementation until explicit confirmation is granted.**

## 5. Exemptions

This protocol does **NOT** apply to:
- Pure investigatory questions ("explain X", "where is Y", "why did Z happen").
- Trivial fixes where the user provides the exact change ("fix this typo", "add this import", "rename this variable").
- Follow-up implementation steps on an already-confirmed understanding lock or approved implementation plan.
- Explicit user override ("just do it", "skip the questions", "I know what I want").

## Reference

- **See skill:** [grill-me](file:///home/isurwars/Projects/Correlation/.agents/skills/grill-me/SKILL.md) for the full interrogation workflow, decision log template, and exit criteria.
