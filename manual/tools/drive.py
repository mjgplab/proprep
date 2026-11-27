#!/usr/bin/env python3
"""Drive an interactive terminal program (ProPrep) from a non-interactive shell.

Each session runs the program inside a pseudo-terminal, owned by a small
background daemon. Commands talk to it through files in the session folder:

  drive.py start NAME [--cwd DIR] [--wait S] [-- CMD ARGS...]
        Launch CMD (default: ProPrep env `proprep`) in a PTY; print first screen.
  drive.py send NAME TEXT [--wait S] [--timeout S] [--no-enter]
        Type TEXT (Enter is added unless --no-enter), then print everything the
        program printed until it has been quiet for --wait seconds (default 2)
        or --timeout (default 120) passes.
        Escapes in TEXT: \\e (ESC) \\t \\x03 (Ctrl-C) \\x04 (Ctrl-D), and
        the words <UP> <DOWN> <LEFT> <RIGHT> <ENTER> <CTRL-C>.
  drive.py read NAME [--wait S] [--timeout S]
        Print new output without typing anything (for long-running steps).
  drive.py screen NAME [--lines N]
        Print the last N lines (default 60) of the whole transcript.
  drive.py status NAME      Is the program still running?
  drive.py stop NAME        Kill the program and the daemon.

Transcripts: SESSION/transcript.txt (ANSI stripped, includes '>>> SENT:' marks)
and SESSION/raw.log. Sessions live under manual/runs/_sessions/NAME/.
"""
import argparse
import fcntl
import os
import pty
import re
import select
import signal
import struct
import sys
import termios
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "runs", "_sessions"))
DEFAULT_CMD = ["/opt/homebrew/Caskroom/miniforge/base/envs/ProPrep/bin/proprep"]
ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]|\x1b\][^\x07]*(\x07|\x1b\\)|\x1b[()][A-Za-z0-9]|\x1b[=>]")


def sdir(name):
    return os.path.join(ROOT, name)


def strip(b: bytes) -> str:
    t = b.decode("utf-8", "replace")
    t = ANSI.sub("", t)
    t = t.replace("\r\n", "\n").replace("\r", "\n")
    return re.sub(r"[\x00-\x08\x0b-\x1f\x7f]", "", t)


def daemon(name, cmd, cwd):
    d = sdir(name)
    pid, fd = pty.fork()
    if pid == 0:
        os.chdir(cwd)
        env = dict(os.environ)
        env.update({"TERM": "xterm-256color", "COLUMNS": "140", "LINES": "50",
                    "PYTHONUNBUFFERED": "1"})
        os.execvpe(cmd[0], cmd, env)
    fcntl.ioctl(fd, termios.TIOCSWINSZ, struct.pack("HHHH", 50, 140, 0, 0))
    with open(os.path.join(d, "child.pid"), "w") as f:
        f.write(str(pid))
    inq = os.path.join(d, "inbox")
    raw = open(os.path.join(d, "raw.log"), "ab", buffering=0)
    txt = open(os.path.join(d, "transcript.txt"), "a", buffering=1)
    alive = True
    while alive:
        r, _, _ = select.select([fd], [], [], 0.1)
        if fd in r:
            try:
                data = os.read(fd, 65536)
            except OSError:
                data = b""
            if not data:
                alive = False
            else:
                raw.write(data)
                txt.write(strip(data))
        # pending input files: inbox/NNNN.in
        try:
            items = sorted(os.listdir(inq))
        except FileNotFoundError:
            items = []
        for it in items:
            p = os.path.join(inq, it)
            with open(p, "rb") as f:
                payload = f.read()
            os.remove(p)
            txt.write(f"\n>>> SENT: {payload!r}\n")
            try:
                os.write(fd, payload)
            except OSError:
                alive = False
        wp = os.waitpid(pid, os.WNOHANG)
        if wp[0] != 0:
            # drain
            try:
                while True:
                    data = os.read(fd, 65536)
                    if not data:
                        break
                    raw.write(data)
                    txt.write(strip(data))
            except OSError:
                pass
            alive = False
            with open(os.path.join(d, "exit_status"), "w") as f:
                f.write(str(wp[1] >> 8))
    txt.write("\n>>> PROGRAM EXITED\n")


def tsize(name):
    p = os.path.join(sdir(name), "transcript.txt")
    return os.path.getsize(p) if os.path.exists(p) else 0


def wait_quiet(name, start, wait, timeout):
    p = os.path.join(sdir(name), "transcript.txt")
    t0 = time.time()
    last = tsize(name)
    last_change = time.time()
    while True:
        time.sleep(0.25)
        s = tsize(name)
        if s != last:
            last, last_change = s, time.time()
        exited = os.path.exists(os.path.join(sdir(name), "exit_status"))
        if (time.time() - last_change >= wait and s > start) or exited:
            if exited:
                time.sleep(0.5)
            break
        if time.time() - t0 > timeout:
            print(f"[drive] timeout after {timeout}s (program may still be working; use `read`)")
            break
    with open(p, "rb") as f:
        f.seek(start)
        out = f.read().decode("utf-8", "replace")
    # collapse runs of blank lines
    out = re.sub(r"\n{3,}", "\n\n", out)
    print(out)
    with open(os.path.join(sdir(name), "read.mark"), "w") as f:
        f.write(str(tsize(name)))
    if os.path.exists(os.path.join(sdir(name), "exit_status")):
        print("[drive] program has exited")


def decode_text(text, enter):
    for k, v in {"<UP>": "\x1b[A", "<DOWN>": "\x1b[B", "<RIGHT>": "\x1b[C", "<LEFT>": "\x1b[D",
                 "<ENTER>": "\r", "<CTRL-C>": "\x03", "<CTRL-D>": "\x04", "<ESC>": "\x1b"}.items():
        text = text.replace(k, v)
    text = text.replace("\\e", "\x1b").replace("\\t", "\t").replace("\\x03", "\x03").replace("\\x04", "\x04")
    if enter:
        text += "\r"
    return text.encode()


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="op", required=True)
    s = sub.add_parser("start"); s.add_argument("name"); s.add_argument("--cwd")
    s.add_argument("--wait", type=float, default=4); s.add_argument("--timeout", type=float, default=120)
    s.add_argument("cmd", nargs="*")
    s = sub.add_parser("send"); s.add_argument("name"); s.add_argument("text")
    s.add_argument("--wait", type=float, default=2); s.add_argument("--timeout", type=float, default=120)
    s.add_argument("--no-enter", action="store_true")
    s = sub.add_parser("read"); s.add_argument("name")
    s.add_argument("--wait", type=float, default=3); s.add_argument("--timeout", type=float, default=300)
    s = sub.add_parser("screen"); s.add_argument("name"); s.add_argument("--lines", type=int, default=60)
    s = sub.add_parser("status"); s.add_argument("name")
    s = sub.add_parser("stop"); s.add_argument("name")
    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        i = argv.index("--"); extra = argv[i+1:]; argv = argv[:i]
    a = ap.parse_args(argv)
    if a.op == "start":
        a.cmd = extra
    d = sdir(a.name)

    if a.op == "start":
        if os.path.exists(os.path.join(d, "child.pid")) and not os.path.exists(os.path.join(d, "exit_status")):
            sys.exit(f"session {a.name} already running; stop it first")
        os.makedirs(os.path.join(d, "inbox"), exist_ok=True)
        for f in ("exit_status", "child.pid"):
            if os.path.exists(os.path.join(d, f)):
                os.remove(os.path.join(d, f))
        cmd = [c for c in a.cmd if c != "--"] or DEFAULT_CMD
        cwd = os.path.abspath(a.cwd or os.getcwd())
        os.makedirs(cwd, exist_ok=True)
        start = tsize(a.name)
        with open(os.path.join(d, "transcript.txt"), "a") as f:
            f.write(f"\n>>> START: {cmd} in {cwd} at {time.ctime()}\n")
        start = tsize(a.name)
        if os.fork() == 0:
            os.setsid()
            if os.fork() == 0:
                devnull = os.open(os.devnull, os.O_RDWR)
                for x in (0, 1, 2):
                    os.dup2(devnull, x)
                with open(os.path.join(d, "daemon.pid"), "w") as f:
                    f.write(str(os.getpid()))
                daemon(a.name, cmd, cwd)
            os._exit(0)
        wait_quiet(a.name, start, a.wait, a.timeout)
    elif a.op == "send":
        if os.path.exists(os.path.join(d, "exit_status")):
            sys.exit("[drive] program has exited; start a new session")
        start = tsize(a.name)
        p = os.path.join(d, "inbox", f"{time.time_ns()}.in")
        with open(p + ".tmp", "wb") as f:
            f.write(decode_text(a.text, not a.no_enter))
        os.rename(p + ".tmp", p)
        wait_quiet(a.name, start, a.wait, a.timeout)
    elif a.op == "read":
        mark = os.path.join(d, "read.mark")
        start = int(open(mark).read()) if os.path.exists(mark) else 0
        wait_quiet(a.name, start, a.wait, a.timeout)
    elif a.op == "screen":
        with open(os.path.join(d, "transcript.txt"), errors="replace") as f:
            print("".join(f.readlines()[-a.lines:]))
    elif a.op == "status":
        print("exited" if os.path.exists(os.path.join(d, "exit_status")) else "running")
    elif a.op == "stop":
        for f in ("child.pid", "daemon.pid"):
            p = os.path.join(d, f)
            if os.path.exists(p):
                try:
                    os.kill(int(open(p).read()), signal.SIGKILL)
                except (ProcessLookupError, ValueError):
                    pass
        open(os.path.join(d, "exit_status"), "w").write("killed")
        print("stopped")


if __name__ == "__main__":
    main()
