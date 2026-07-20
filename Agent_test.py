"""
Local test harness for the "Multi-Omics Thesis Data Engineer" agent.

Runs entirely through the Claude API (no Managed Agents console needed).
Uses the built-in code execution tool so Claude actually runs Python/bash
in a sandbox, rather than a hand-defined function-calling tool.

Because the system prompt requires "one micro-task at a time, then halt,
report shapes/metrics, and wait for explicit written approval," this
script surfaces every tool call to you and asks for confirmation before
sending the result back to Claude. You are the approval gate.
"""

import os
from dotenv import load_dotenv
import anthropic

# Loads variables from a .env file in the same folder as this script,
# regardless of what your shell/IDE kernel does or doesn't inherit.
load_dotenv()

api_key = os.environ.get("ANTHROPIC_API_KEY")
if not api_key:
    raise RuntimeError(
        "ANTHROPIC_API_KEY not found. Make sure you have a .env file in the "
        "same folder as this script containing:\nANTHROPIC_API_KEY=your-key-here"
    )

client = anthropic.Anthropic(api_key=api_key)

MODEL = "claude-opus-4-8"

SYSTEM_PROMPT = """You are a computational biology assistant specializing in single-cell and bulk multi-omics (RNA-seq, ATAC-seq) data engineering, using Scanpy, Pandas, and AnnData (.h5ad) objects. Never use deprecated or hallucinated functions (e.g. no sc.tl.enrichment; use Squidpy sq.gr.nhood_enrichment for spatial/network metrics). You are strictly forbidden from running multi-step pipelines unattended: complete only one micro-task at a time, then halt, report shapes/metrics, and wait for explicit written approval before continuing. Before each code execution, state the packages used, your assumptions about data/biology, and the specific processing target. Begin any new task with a context-inspection audit: environment/library versions, metadata columns (age, sex, race, lifestyle, batch), and unit consistency. Standardize schema to adata.obs['sample'] and adata.obs['cell_type']. Reconcile gene synonyms across datasets. Ground analyses in healthy-tissue biology; never apply cancer-genetics methods (e.g. CNV inference) to healthy samples. Assertive Vetting: When evaluating any new or publicly available datasets for inclusion, you must write explicit, deterministic Boolean tests (e.g., Python assert statements verifying column existence, non-null status, correct data types, and numeric ranges) to validate data integrity before integrating them."""

# code_execution_20260120 = current version, no beta header required,
# supports Opus 4.6+/4.8, Sonnet 4.5+/4.6/5. Also gives you Python REPL
# state persistence across turns within the same conversation.
TOOLS = [
    {"type": "code_execution_20260120", "name": "code_execution"}
]

# The Files API itself still needs this beta header, even though
# code_execution_20260120 doesn't. Any message that references an
# uploaded file (via container_upload) must send it too.
FILES_BETA_HEADER = "files-api-2025-04-14"


def upload_file(path):
    """
    Upload a local file to Anthropic's Files API and return its file_id.
    This does NOT put the file in the sandbox yet — it just stores it on
    Anthropic's side so you can reference it (cheaply, repeatedly) in
    later messages via a container_upload block.
    """
    with open(path, "rb") as f:
        filename = os.path.basename(path)
        uploaded = client.beta.files.upload(
            file=(filename, f),
            extra_headers={"anthropic-beta": FILES_BETA_HEADER},
        )
    print(f"Uploaded {filename} -> file_id={uploaded.id}")
    return uploaded.id


def run_step(messages, container_id=None):
    """
    Send one turn to Claude, print any text/tool activity, return
    (response, container_id). Passing container_id reuses the same
    sandbox container as a prior turn, so uploaded files and any
    Python/bash state survive between messages instead of vanishing
    into a fresh, empty container each time.
    """
    kwargs = dict(
        model=MODEL,
        max_tokens=4096,
        system=SYSTEM_PROMPT,
        tools=TOOLS,
        messages=messages,
        extra_headers={"anthropic-beta": FILES_BETA_HEADER},
    )
    if container_id:
        kwargs["container"] = container_id

    response = client.messages.create(**kwargs)

    for block in response.content:
        if block.type == "text":
            print("\n--- Claude ---")
            print(block.text)
        elif block.type == "server_tool_use":
            print("\n--- Claude wants to run code ---")
            print(block.input.get("code", block.input))
        elif block.type == "code_execution_tool_result":
            print("\n--- Sandbox output ---")
            print(block.content)

    new_container_id = getattr(response, "container", None)
    new_container_id = new_container_id.id if new_container_id else container_id

    return response, new_container_id


def main():
    # -----------------------------------------------------------------
    # OPTIONAL: attach real files here. List local paths to your
    # Morandini 2024 (GSE193140/141/142) files. Leave the list empty to
    # just run the environment audit with no data attached yet.
    #
    # Example:
    # file_paths = [
    #     "/Users/marikaclark/dark_matter_atx/data/GSE193140_counts.csv",
    #     "/Users/marikaclark/dark_matter_atx/data/GSE193141_atac.mtx",
    # ]
    # -----------------------------------------------------------------
    file_paths = []

    content = []
    for path in file_paths:
        file_id = upload_file(path)
        content.append({"type": "container_upload", "file_id": file_id})

    task_text = (
        "Begin the context-inspection audit for my thesis pipeline. "
        + (
            "I've attached the Morandini 2024 ATAC-Clock dataset "
            "(GSE193140/141/142) files. Inspect what's actually in them "
            "(file type, shape, columns/headers) before assuming structure. "
            if file_paths
            else "I'll upload the Morandini 2024 ATAC-Clock dataset "
            "(GSE193140/141/142) files separately. For now, just check "
            "the environment: confirm scanpy, pandas, and anndata "
            "versions available in this sandbox. "
        )
    )
    content.append({"type": "text", "text": task_text})

    messages = [{"role": "user", "content": content}]

    response, container_id = run_step(messages)
    messages.append({"role": "assistant", "content": response.content})

    # Ongoing conversation loop. After every Claude turn, you type your
    # actual reply (approval, an answer to a question, a correction,
    # whatever) and it gets sent back as the next user message. This is
    # your approval gate: Claude only proceeds when you type something.
    # container_id is threaded through so the sandbox (and any files
    # uploaded into it) persists across turns instead of expiring.
    print("\n" + "=" * 60)
    while True:
        user_reply = input("\nYour reply (or 'quit' to stop): ").strip()

        if user_reply.lower() in ("quit", "exit"):
            print("Stopped.")
            break

        messages.append({"role": "user", "content": user_reply})
        response, container_id = run_step(messages, container_id=container_id)
        messages.append({"role": "assistant", "content": response.content})
        print("\n" + "=" * 60)


if __name__ == "__main__":
    main()