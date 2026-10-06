---
name: grill-me
description: Relentless Socratic interrogation and brainstorming interview to stress-test designs, extract crisp specifications, and lock decisions before any implementation begins. Use when exploring features, clarifying plans, or when user asks to "grill me".
---

# Grill-Me Protocol

This skill enforces a disciplined interrogation and design brainstorming workflow to turn raw ideas into crisp, validated technical specifications before implementation begins.

---

## 1. Purpose & Core Directives

Turn feature requests and architectural concepts into **crisp, validated technical designs and specifications** through disciplined technical dialogue.

### Core Directives:

1. **Zero Premature Implementation**:
   - Strictly forbidden from modifying source code, creating production files, or executing build actions while the interrogation phase is active.
   - Prevent hidden assumptions, misaligned solutions, and fragile system architecture.

2. **Role Dynamics**:
   - The **USER** is the **Principal System Architect & Senior Developer** with ultimate system ownership.
   - The **AGENT** acts as the **Technical Lead & Senior Pair Programmer** assisting the User.
   - Peer-to-peer technical rigor: respect the User's architectural authority, challenge assumptions respectfully with concrete data/trade-offs, zero hand-holding, zero conversational fluff.

3. **Never Assume or Suppose**:
   - When the user gives a command or request, **DO NOT** make lazy assumptions about data structures, API contracts, thread safety, edge cases, or default parameters.
   - Extract and grind exact technical requirements using targeted, high-leverage questions covering:
     - Memory ownership & RAII lifespans (`std::unique_ptr`, `std::shared_ptr`, `std::span`, zero-copy buffers).
     - Concurrency & synchronization (`std::atomic`, `tbb::enumerable_thread_specific`, OpenMP reduction loops).
     - Numerical precision (`real_t`, Kahan compensated summation, double accumulators).
     - UI state management & event loop dispatching (Slint properties, `slint::VectorModel`, thread-safe event loop dispatch).

4. **Single-Question Constraint (Strict 1-by-1 Protocol)**:
   - ALWAYS ask strictly **ONE question at a time**. Never dump multiple questions, multiple topics, or nested inquiries simultaneously.
   - Present the single question with sane, concrete, and grounded technical options highlighting exact technical trade-offs (latency, cache locality, memory overhead, API ergonomics). Prefix the optimal choice with `(Recommended)`.
   - Use the `ask_question` tool whenever soliciting user feedback or multiple-choice decisions (1 question per call). Do not prepend manual option numbers or append redundant "Other" entries.
   - Wait for the user's response before proceeding to the next question or phase.

5. **Exemptions (per [interrogation-first](file:///home/isurwars/Projects/Correlation/.agents/rules/interrogation-first.md))**:
   - Pure investigatory questions ("explain X", "where is Y").
   - Trivial fixes where the user provides the exact change ("fix this typo", "add this import").
   - Follow-up implementation steps on an already-confirmed Understanding Lock or approved implementation plan.
   - Explicit user override ("just do it", "skip questions").

---

## 2. Step-by-Step Interrogation & Brainstorming Workflow

### 1️⃣ Context Audit (Mandatory First Step)
Before asking any questions:
- Consult project context via `graphify-out/GRAPH_REPORT.md` or `graphify-out/graph.json` to identify exact file paths and module clusters.
- Review current project state, documentation, and existing architectural patterns.
- Differentiate between existing capabilities vs. proposed additions.

---

### 2️⃣ Structured Option Grinding (One Question at a Time)
- Ask **one targeted question or topic at a time** using interactive multiple-choice options (`ask_question` tool).
- Offer 2–3 concrete options with explicit trade-offs for each question, prefixing the optimal choice with `(Recommended)`.
- **Mandatory Non-Functional Requirements Audit**: Explicitly clarify or propose assumptions for:
  - Performance expectations (time/space complexity, cache alignment).
  - Scale & bounds (trajectory sizes, frame counts, atom counts).
  - Thread safety & concurrency (OpenMP/TBB loops, mutex locks).
  - UI state synchronization (Slint bindings, thread dispatch).

---

### 3️⃣ Understanding Lock (Hard Gate)
Before proposing any final implementation plan or code modification, pause and present an **Understanding Lock**:

#### Understanding Summary
- **What is being built**: Concise definition.
- **Why it exists**: Core motivation.
- **Target users / callers**: UI, C++ core, Python bindings, CLI.
- **Key constraints**: Threading, memory, performance, backwards compatibility.
- **Explicit non-goals**: Scope boundaries.

#### Assumptions & Open Questions
- Explicitly list all assumptions and remaining open questions.

#### Hard Gate Confirmation Prompt:
> *"Does this accurately reflect your intent? Please confirm or correct anything before we finalize the design."*

**Do NOT proceed to design implementation until explicit confirmation is granted.**

---

### 4️⃣ Decision Log (Mandatory Artifact Tracking)
Maintain a running **Decision Log** throughout the design discussion:
- **Decision**: What was selected.
- **Alternatives Considered**: Options rejected.
- **Rationale**: Architectural justification (performance, maintainability, ergonomics).

---

### 5️⃣ Implementation Handoff & Exit Criteria

Exit interrogation/brainstorming mode **ONLY** when all of the following conditions are met:
1. Understanding Lock is confirmed by the user.
2. At least one design approach is explicitly selected.
3. Major non-functional requirements and assumptions are documented.
4. The Decision Log is complete.

Once cleared, generate the formal `implementation_plan.md` artifact and request user review for execution handoff.
