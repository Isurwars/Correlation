---
name: handoff
description: Compress and document conversation state, active technical decisions, and pending tasks into a structured handoff artifact. Use when ending a session, handing work off to another agent/session, or saving progress.
---

# Session Handoff

Produce a dense, high-signal handoff summary to preserve momentum, architectural decisions, and next immediate steps across sessions or agents.

## Phase 1: Context Distillation

Collect the core vectors of the active session:
1. **Goal Achieved:** What was requested and what was completed.
2. **Current State:** Git status, passing/failing tests, build status.
3. **Key Decisions Made:** Architectural or algorithmic choices locked with the user.
4. **Active Blockers / Gaps:** Unresolved questions or technical obstacles.
5. **Next Immediate Steps:** The precise, atomic next actions for the incoming agent/session.

## Phase 2: Handoff Template

Output the summary using the standard format:

```markdown
# Session Handoff Summary

## 1. Executive Summary
- **Current Objective:** [Brief statement of overarching task]
- **Status:** [IN_PROGRESS / BLOCKED / READY_FOR_VERIFICATION]
- **Active Git Commit/Branch:** [Branch name or HEAD short SHA]

## 2. Work Completed
- [Completed item 1 with file links]
- [Completed item 2 with file links]

## 3. Architecture & Decision Log
| Decision | Rationale | Alternatives Rejected |
| :--- | :--- | :--- |
| [Decision 1] | [Rationale] | [Rejected] |

## 4. Immediate Next Actions (Tracer Bullets)
1. **[Next Step 1]:** [File and exact action]
2. **[Next Step 2]:** [Test or build command to execute]

## 5. Critical Invariants & Warnings
- [Known gotchas, thread constraints, or sensitive files]
```
