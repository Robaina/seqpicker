"""Utility functions and classes for the reducefasta package."""

from __future__ import annotations

import random
import shutil
import string
import subprocess
from pathlib import Path
from typing import Optional


class TemporaryFilePath:
    """Context manager for creating temporary files that are automatically cleaned up."""

    def __init__(
        self,
        work_dir: Optional[Path] = None,
        extension: Optional[str] = None,
    ):
        """Initialize temporary file path manager.

        Args:
            work_dir: Working directory for the temporary file
            extension: Optional file extension
        """
        self.work_dir = Path(work_dir).resolve() if work_dir else Path.cwd()
        self.extension = extension or ""

    def __enter__(self) -> Path:
        temp_id = "".join(random.choices(string.ascii_lowercase, k=10))
        temp_file_name = f"temp_{temp_id}{self.extension}"
        self.file_path = self.work_dir / temp_file_name
        return self.file_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file_path.unlink(missing_ok=True)


class TemporaryDirectoryPath:
    """Context manager for creating temporary directories that are automatically cleaned up."""

    def __init__(self, work_dir: Optional[Path] = None):
        """Initialize temporary directory manager.

        Args:
            work_dir: Parent directory for the temporary directory
        """
        self.work_dir = Path(work_dir).resolve() if work_dir else Path.cwd()

    def __enter__(self) -> Path:
        temp_id = "".join(random.choices(string.ascii_lowercase, k=10))
        self.dir_path = self.work_dir / f"temp_{temp_id}"
        self.dir_path.mkdir(parents=True, exist_ok=True)
        return self.dir_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.dir_path.exists():
            shutil.rmtree(self.dir_path)


def run_command(
    cmd: str,
    suppress_output: bool = False,
    work_dir: Optional[Path] = None,
) -> subprocess.CompletedProcess:
    """Execute a command in the shell.

    Args:
        cmd: Command string to execute
        suppress_output: Whether to suppress command output
        work_dir: Working directory for command execution

    Returns:
        CompletedProcess instance with command execution results
    """
    if suppress_output:
        cmd = f"{cmd} >/dev/null 2>&1"

    return subprocess.run(
        cmd,
        shell=True,
        cwd=work_dir,
        check=True,
        text=True,
    )
