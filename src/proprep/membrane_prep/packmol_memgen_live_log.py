"""
Start packmol-memgen with its log files written line by line.

Run as a script, by the interpreter of packmol-memgen's own shebang line:

    <python> packmol_memgen_live_log.py /path/to/packmol-memgen <its arguments>

packmol-memgen copies packmol's output into packmol.log through an ordinary
file object, which Python flushes 8 KB at a time. One packing loop writes about
1 KB, so the log, and with it ProPrep's progress display, moved once every five
or six loops: every five minutes on a system whose loops take 40 s. Measured on
a 70 A POPC patch: the log grew 2 times in 40 s as it was, 50 times with this.

Nothing about packmol-memgen is changed except the buffering of text files it
opens for writing whose name ends in ".log". It imports nothing from ProPrep,
because the interpreter that runs packmol-memgen need not have ProPrep.
"""

import builtins
import runpy
import sys

_open = builtins.open


def _line_buffered_logs(file, mode="r", buffering=-1, *args, **kwargs):
    if (isinstance(file, str) and file.endswith(".log") and "w" in mode
            and "b" not in mode and buffering == -1):
        buffering = 1
    return _open(file, mode, buffering, *args, **kwargs)


if __name__ == "__main__":
    builtins.open = _line_buffered_logs
    executable = sys.argv[1]
    sys.argv = sys.argv[1:]
    runpy.run_path(executable, run_name="__main__")
