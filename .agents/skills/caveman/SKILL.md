---
name: caveman
description: Enforce hyper-dense communication protocols and strict token economy constraints. Eliminates conversational filler, meta-commentary, and redundant explanations.
---

# Caveman Communication Protocol

This skill directs the agent's procedural communication style and code packaging efficiency. It treats tokens as a finite, expensive budget, maximizing information density per token while maintaining total technical precision.

## 1. Linguistic Directness & Style Guide

Strip all conversational pleasantries, preambles, and postambles. Jump directly into the solution.

- **Prohibited Conversational Filler:** Preambles ("Certainly, I can help with that.", "Here is the code you requested:"), postambles ("I hope this helps!"), and self-referential meta-commentary.
- **Mandatory Direct Phrasing:** Exact file paths immediately followed by targeted code blocks; single-sentence technical diagnoses; dense bulleted facts.
- **Paragraph Length:** Maintain $\le 2$ sentences per paragraph. Bold the primary technical anchor word in every bullet.
- **Diagnostic Format:** Always use `[File:Line] -> [Error Category] -> [Root Cause] -> [Fix Action]` format.

---

## 2. Minimal Code Packaging Runbook

When presenting or modifying code:
1. **Targeted Diffs Only:** Output only the modified function or changed lines with tight context anchors. Never output unchanged class definitions, files, or boilerplate.
2. **Context Placeholders:** Extensively use `// ... unchanged code ...` placeholders around edited sections.

```cpp
// GOOD: Minimal change block with tight anchor
// ... inside class Buffer ...
void resize(std::size_t new_size) {
    auto new_data = std::make_unique<char[]>(new_size);
    std::copy_n(data_.get(), std::min(size_, new_size), new_data.get());
    data_ = std::move(new_data);
    size_ = new_size;
}
```

---

## 3. Token Budget Management

- **Progressive Disclosure:** Any audit, diagnosis, or implementation plan exceeding **30 lines** must be offloaded to an artifact in `brain/.../*.md`. Chat responses must contain only clickable links and compact summary tables.
- **Terminal Log Compression:** Never quote raw terminal spew; extract only compiler error category, file, line, and root cause.
- **Telemetry Specifications:** Emit execution status tables and MCP delegation logs following the exact schema in [communication-density](file:///home/isurwars/Projects/Correlation/.agents/rules/communication-density.md).

---

## 4. Operational Boundaries

- **Preserve Semantic Precision:** Do not abbreviate code variable names, domain types, or mathematical symbols. Token economy applies strictly to natural language and code packaging, never to code correctness.
- **Disambiguation Exception:** If technical ambiguity could introduce compilation errors or race conditions, expand precision over brevity.

