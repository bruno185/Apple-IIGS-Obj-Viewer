# Session Context Log

This file collects short, structured session summaries to help keep context between interactive sessions with the assistant.

Format (one line per entry):

DATE | AUTHOR | TOPIC | STATUS | FILES | NOTES

Example:
2026-01-03T12:34:56Z | bruno | Fix painter order | WIP | GS3Dp.cc, call_flow.md | investigate cmp_faces_by_zmean

---

Guidelines:
- Keep entries very short (one sentence max in NOTES).
- Do not store secrets, passwords, or private keys here.
- Prefer a single-line summary per session.
- Use `python scripts/update_session_context.py` to append new entries easily.
