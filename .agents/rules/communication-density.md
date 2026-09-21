# Rule: Communication Density & Token Economy

*Activation Mode: Universal / Default Baseline*

## 1. Zero Conversational Filler

- No preambles ("Certainly, I can help with that.").
- No postambles ("I hope this helps! Let me know if you have questions.").
- No meta-commentary about what the agent is about to do or just did.
- Jump directly into the solution, diagnosis, or question.

## 2. Minimal Code Output

- Use `// ... unchanged code ...` placeholders extensively. Never print untouched structural logic or boilerplate code blocks.
- When modifying existing code, output only the target function or modified lines with tight context anchors — not the entire class or file.

## 3. Strict Diagnostic Format

When reporting compilation issues, failures, or code review findings, use:

```
[File:Line] → [Error Type] → [Fix Action]
```

**Example:** `[src/buffer.cpp:14] → Linker Error (Unresolved external symbol) → Add target_link_libraries in CMakeLists.txt.`

## 4. Structural Scannability

- Paragraphs must be **≤ 2 sentences**.
- Prefer **tables** for multi-variable comparisons.
- **Bold** the primary technical anchor word in every bullet point.

## 5. Token Economy

| Resource Type      | Waste Vector                   | Mitigation                                                                |
| :----------------- | :----------------------------- | :------------------------------------------------------------------------ |
| **Output Tokens**  | Code repetition / Explanations | Never repeat unchanged code. Use `// ...` placeholders.                   |
| **Input Tokens**   | Bloated error logs             | Extract only file, line, and diagnostic message from terminal output.     |
| **Context Window** | Outdated design drafts         | Overwrite old drafts when a design reaches finality.                      |

## Reference

- **See skill:** [caveman](file:///home/isurwars/Projects/Correlation/.agents/skills/caveman/SKILL.md) for extended communication protocol, code packaging rules, and limitations.
