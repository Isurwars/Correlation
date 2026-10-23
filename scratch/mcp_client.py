import json
import os
import sys
import urllib.request

LLAMA_SERVER_URL = os.environ.get("LLAMA_SERVER_URL", "http://127.0.0.1:9931")

def ask_model(prompt: str, system_prompt: str = "You are a specialized C++23 execution worker.") -> str:
    url = f"{LLAMA_SERVER_URL.rstrip('/')}/v1/chat/completions"
    payload = {
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": prompt}
        ],
        "temperature": 0.2,
        "max_tokens": 4096
    }
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30.0) as resp:
            res = json.loads(resp.read().decode('utf-8'))
            msg = res["choices"][0]["message"]
            content = msg.get("content", "")
            return content
    except Exception as e:
        sys.stderr.write(f"Error querying llama-server: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'r') as f:
            prompt_text = f.read()
    else:
        prompt_text = sys.stdin.read()
    print(ask_model(prompt_text))
