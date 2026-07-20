---
name: caveman-navigation
description: Minimizes token overhead during architectural analysis and codebase queries.
---

## Core Operational Constraints
1. **Zero Exploratory Reading:** You are forbidden from opening a source file unless its path is explicitly mentioned in `graphify-out/GRAPH_REPORT.md` as a structural dependency of the feature being modified.
2. **Context Budgeting:** Treat every 1,000 tokens like physical currency. Do not pull in header files or module files unless an explicit compile error or link dependency demands it.
3. **Graph Leverage:** Use `graphify-out/graph.json` or the markdown report to answer questions about "What modules talk to each other?" never grep code imports to figure this out.

## Execution Pattern
- **Step 1:** Open `graphify-out/GRAPH_REPORT.md`. Identify the "Community Cluster" or target node.
- **Step 2:** Formulate a single, direct file path to inspect.
- **Step 3:** Perform the modification or analysis on that single node, then immediately flush non-essential context buffers.