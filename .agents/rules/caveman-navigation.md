---
trigger: always_on
---

# Rule: Caveman Navigation & Architecture Mapping
*Activation Mode: Always On*

## Core Constraint (Context Diet)
- **Prohibited Action:** Opening source code files directly as a baseline discovery action is strictly prohibited. Do not conduct broad file searches or recursive grep queries across the workspace blindly.
- **First-Action Enforce:** Every new conversational turn, distinct sub-task, or multi-file architectural analysis MUST begin by reading the structural layout database artifact located at:
  `graphify-out/GRAPH_REPORT.md`

## Tool Sequence
1. Consult `graphify-out/GRAPH_REPORT.md` and `graphify-out/graph.json` first to look up dependency links, file clusters, and component definitions.
2. Target specific lines or isolated code components *only* after locating their precise coordinates in the graph report.