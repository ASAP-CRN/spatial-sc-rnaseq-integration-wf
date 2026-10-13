import sys
from datetime import datetime
from pathlib import Path


class _Tee:
    """Write to both a file and the original stream simultaneously."""

    def __init__(self, stream, log_fh):
        self._stream = stream
        self._log = log_fh

    def write(self, data):
        self._stream.write(data)
        self._log.write(data)

    def flush(self):
        self._stream.flush()
        self._log.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)


def log_run_command(log_file: Path) -> None:
    """
    tee calling script invocation with a timestamp to `log_file`.

    Parameters
    ----------
    log_file : Path
        Log file path. Created (with parent dirs) if it does not exist.

    Returns
    -------
    None
    """
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    fh = open(log_file, "w", buffering=1)  # line-buffered

    cmd = " ".join(sys.argv)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fh.write(f"[{timestamp}] {cmd}\n")

    sys.stdout = _Tee(sys.__stdout__, fh)
    sys.stderr = _Tee(sys.__stderr__, fh)
