#!/usr/bin/env python3
"""
MCP Worker client script for Correlation agent.
Dispatches scoped prompts to the local llama-server following the mcp-delegation protocol.
"""

import sys
import json
import urllib.request
import urllib.error

LLAMA_URL = "http://127.0.0.1:9931/v1/chat/completions"

def query_mcp(prompt: str, system_prompt: str = "You are a specialized C++23 execution worker. Output clean, optimal C++23 code only without conversational preamble.", timeout: float = 45.0) -> str:
    payload = {
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": prompt}
        ],
        "temperature": 0.2,
        "top_p": 0.95,
        "max_tokens": 1024
    }
    
    req = urllib.request.Request(
        LLAMA_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            res_data = json.loads(response.read().decode("utf-8"))
            choice = res_data.get("choices", [{}])[0]
            msg = choice.get("message", {})
            content = msg.get("content", "").strip()
            if not content and "reasoning_content" in msg:
                content = msg["reasoning_content"].strip()
            return content
    except Exception as e:
        print(f"MCP_ERROR: {e}", file=sys.stderr)
        return ""

import argparse

def main():
    parser = argparse.ArgumentParser(description="MCP worker client")
    parser.add_argument("--prompt", type=str, default="", help="Prompt string")
    parser.add_argument("--system", type=str, default="You are a specialized C++23 execution worker. Output clean, optimal C++23 code only without conversational preamble.", help="System prompt")
    parser.add_argument("--file", type=str, default="", help="Path to prompt file")
    parser.add_argument("--timeout", type=float, default=60.0, help="Request timeout in seconds")
    args = parser.parse_args()

    user_prompt = args.prompt
    if args.file:
        with open(args.file, "r", encoding="utf-8") as f:
            user_prompt = f.read()
    elif not user_prompt and not sys.stdin.isatty():
        user_prompt = sys.stdin.read()

    if not user_prompt:
        print("Error: No prompt provided", file=sys.stderr)
        sys.exit(1)

    result = query_mcp(user_prompt, system_prompt=args.system, timeout=args.timeout)
    if result:
        print(result)
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
