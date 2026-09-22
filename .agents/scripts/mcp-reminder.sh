#!/usr/bin/env bash
set -euo pipefail

# Read stdin into variable to ensure valid JSON input consumption
cat > /dev/null

cat <<'EOF'
{
  "injectSteps": [
    {
      "ephemeralMessage": "[MANDATORY MCP INVARIANT] You MUST invoke `call_mcp_tool` targeting `llama-cpp-mcp` (`ask_model` or `analyse_project`) as your FIRST tool call before answering architectural questions, designing implementations, or modifying files. Direct execution without MCP consultation is strictly prohibited unless `llama-cpp-mcp` returns an error or is unreachable."
    }
  ]
}
EOF
