"""
GC-NKGraph-Atlas common logging utilities.

Centralized logging with NEEDS_REVIEW markers, dataset failure tracking,
and structured log output to results/logs/LOG.md.
"""

import sys
import os
from datetime import datetime
from pathlib import Path
from typing import Optional


LOG_ENTRY_TEMPLATE = """---
## {timestamp} | {phase} | {status}

**File:** `{script}`
**Config:** `{config}`

{message}

"""


class Logger:
    """Project-level logger writing structured markdown to results/logs/LOG.md."""

    def __init__(self, log_path: str = "results/logs/LOG.md", project_root: Optional[str] = None):
        if project_root:
            self.log_path = Path(project_root) / log_path
        else:
            self.log_path = Path(log_path)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        if not self.log_path.exists():
            self.log_path.write_text("# GC-NKGraph-Atlas Execution Log\n\n")

    def _write_entry(self, phase: str, status: str, message: str, script: str = "", config: str = ""):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
        entry = LOG_ENTRY_TEMPLATE.format(
            timestamp=timestamp,
            phase=phase,
            status=status,
            script=script or sys.argv[0],
            config=config or "N/A",
            message=message,
        )
        with open(self.log_path, "a", encoding="utf-8") as f:
            f.write(entry)

    def ok(self, phase: str, message: str, script: str = "", config: str = ""):
        self._write_entry(phase, "OK", message, script, config)

    def skip(self, phase: str, message: str, script: str = "", config: str = ""):
        self._write_entry(phase, "SKIPPED", message, script, config)

    def fail(self, phase: str, message: str, script: str = "", config: str = ""):
        self._write_entry(phase, "FAILED", message, script, config)

    def needs_review(self, phase: str, message: str, script: str = "", config: str = ""):
        """Use for biologically ambiguous decisions requiring human review."""
        self._write_entry(phase, "NEEDS_REVIEW", message, script, config)

    def data_unavailable(self, phase: str, dataset: str, message: str, script: str = ""):
        self._write_entry(phase, "DATA_UNAVAILABLE", f"**Dataset:** {dataset}\n\n{message}", script)

    def add_dataset_failure(self, dataset: str, reason: str):
        """Log a failed dataset with reason to the log."""
        self._write_entry("DATA", "DATASET_FAILED", f"**Dataset:** {dataset}\n**Reason:** {reason}")


class NeedsReviewError(Exception):
    """Raise when a biologically critical decision needs human review."""
    pass


def make_logger(project_root: Optional[str] = None) -> Logger:
    return Logger(project_root=project_root)
