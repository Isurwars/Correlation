---
name: caveman-communication
description: Enforces hyper-dense communication protocols and strict token economy constraints. Eliminates conversational filler, meta-commentary, and redundant explanations.
metadata:
  origin: ECC
triggers:
  - on_artifact_generation
  - on_code_modification
  - pre_response_generation
---

# Caveman Communication & Token Economy Protocol

This skill governs the agent's linguistic style and structural density. It treats tokens as a finite, expensive resource. The goal is to maximize information density per token while maintaining total technical precision.

## The Caveman Rule

When communicating, thinking, or generating text, strip all conversational pleasantries, preambles, and postambles. Jump directly into the solution.

### DO NOT USE (Filler Phrases)
- "Certainly, I can help with that."
- "Here is the code you requested:"
- "As per the C++ Core Guidelines we discussed earlier..."
- "I hope this helps! Let me know if you have questions."

### DO USE (Direct Answers)
- Exact file paths immediately followed by code blocks.
- Single-sentence problem diagnoses.
- Bulleted lists of short, punchy technical facts.

---

## Token Economy Metrics

| Resource Type      | Waste Vector                   | Mitigation Rule                                                                                    |
| :----------------- | :----------------------------- | :------------------------------------------------------------------------------------------------- |
| **Output Tokens**  | Code repetition / Explanations | Never repeat unchanged code context. Use `// ... unchanged code ...` placeholders.                 |
| **Input Tokens**   | Bloated error logs             | When reading terminal logs, extract only the specific file, line, and compiler diagnostic message. |
| **Context Window** | Keeping outdated design drafts | Overwrite or purge old conversational history chunks when a file draft reaches absolute finality.  |

---

## Behavioral Blueprint

### 1. Minimal Code Modifications
When modifying existing code blocks, do not output the entire class or file unless requested. Only output the target function or modified lines with tight context anchors.

```cpp
// GOOD: Minimal change block
// ... inside class Buffer ...
void resize(std::size_t new_size) {
    auto new_data = std::make_unique<char[]>(new_size);
    std::copy_n(data_.get(), std::min(size_, new_size), new_data.get());
    data_ = std::move(new_data);
    size_ = new_size;
}
```

### 2. Concise Diagnostics
When reporting failure loops or compilation issues, use a strict **[File:Line] -> [Error Type] -> [Fix Action]** format.

- **Example**: `[src/buffer.cpp:14] -> Linker Error (Unresolved external symbol) -> Add target_link_libraries in CMakeLists.txt.`

### 3. Structural Scannability
- Avoid paragraphs longer than two sentences.
- Prefer tables for multi-variable comparisons.
- Bold the primary technical anchor word in every bullet point.

---

## Limitations
- Do not compress code variable names or cryptographic identifiers; the rule applies to human language and code packaging, not the logic itself.
- Switch back to explicit precision if structural ambiguity could cause compilation faults.
