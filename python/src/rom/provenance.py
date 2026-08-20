"""Small provenance helpers for generated ROM results."""

from __future__ import annotations

import hashlib
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def _git_value(args: list[str]) -> str:
    try:
        return subprocess.check_output(args, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def file_fingerprint(path: Path, hash_limit_bytes: int = 50_000_000) -> dict:
    path = Path(path)
    if not path.exists():
        return {"path": str(path), "exists": False, "sha256": "unknown"}
    stat = path.stat()
    result = {
        "path": str(path),
        "exists": True,
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if stat.st_size <= hash_limit_bytes:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    else:
        result["sha256"] = "skipped_large_file"
        result["hash_limit_bytes"] = int(hash_limit_bytes)
    return result


def provenance(input_file: Path | None = None, hash_limit_bytes: int = 50_000_000) -> dict:
    status = _git_value(["git", "status", "--porcelain"])
    data = {
        "generation_timestamp": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "git_commit": _git_value(["git", "rev-parse", "HEAD"]),
        "dirty_worktree": bool(status and status != "unknown"),
    }
    if input_file is not None:
        data["input_file"] = file_fingerprint(Path(input_file), hash_limit_bytes)
    return data
