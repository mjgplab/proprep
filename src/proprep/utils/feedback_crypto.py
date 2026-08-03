"""
Encryption helpers for the ProPrep feedback command.

Session transcripts attached to feedback are encrypted to the maintainer's
public key using a libsodium *sealed box* (anonymous X25519 key exchange +
XSalsa20-Poly1305). A sealed box needs only the recipient's public key to
encrypt and only the matching private key to decrypt, with no sender identity
involved. That lets the ciphertext ride inside a *public* GitHub issue while
staying readable exclusively by whoever holds the private key.

PyNaCl is an optional dependency. If it is not installed, or no maintainer
public key has been embedded, encryption is simply unavailable and the caller
omits session context rather than publishing a transcript in the clear.

Maintainer setup (one time):
    proprep --generate-feedback-keypair
paste the printed public key into FEEDBACK_PUBLIC_KEY_B64 below, and keep the
private key it saves to ~/.proprep/feedback_private_key. To read an encrypted
block from an issue later:
    proprep --decrypt-feedback <file-or-'-'>
"""

import base64
import os
import textwrap
from pathlib import Path
from typing import Optional, Tuple

# Maintainer feedback public key (base64-encoded 32-byte X25519 public key).
# Empty string => encrypted session context is disabled (feedback still works;
# the transcript is just never attached). Populate via
# `proprep --generate-feedback-keypair`.
FEEDBACK_PUBLIC_KEY_B64 = "X5r5Jr7t5e9DfUmv9JCuE/dhXpHyckr8FMW+WRBm4Dc="

# Where the maintainer's private key lives (override for testing via env var).
PRIVATE_KEY_PATH = Path(
    os.environ.get("PROPREP_FEEDBACK_KEY",
                   str(Path.home() / ".proprep" / "feedback_private_key"))
)

_ARMOR_HEADER = "-----BEGIN PROPREP ENCRYPTED SESSION-----"
_ARMOR_FOOTER = "-----END PROPREP ENCRYPTED SESSION-----"


def crypto_available() -> bool:
    """True if PyNaCl (the optional crypto dependency) can be imported."""
    try:
        from nacl.public import PrivateKey, PublicKey, SealedBox  # noqa: F401
        return True
    except Exception:
        return False


def public_key_configured() -> bool:
    """True if a maintainer public key has been embedded."""
    return bool(FEEDBACK_PUBLIC_KEY_B64.strip())


def generate_keypair() -> Tuple[str, str]:
    """Generate a fresh X25519 keypair. Returns (public_b64, private_b64)."""
    from nacl.public import PrivateKey
    sk = PrivateKey.generate()
    pub_b64 = base64.b64encode(bytes(sk.public_key)).decode("ascii")
    priv_b64 = base64.b64encode(bytes(sk)).decode("ascii")
    return pub_b64, priv_b64


def encrypt_for_maintainer(plaintext: str,
                           public_key_b64: Optional[str] = None) -> Optional[str]:
    """Seal ``plaintext`` to the maintainer public key.

    Returns an armored, line-wrapped base64 block suitable for pasting into a
    GitHub issue, or ``None`` if crypto is unavailable or no public key is
    configured (so the caller can omit the context instead of leaking it).
    """
    key_b64 = (public_key_b64 if public_key_b64 is not None
               else FEEDBACK_PUBLIC_KEY_B64).strip()
    if not key_b64:
        return None
    try:
        from nacl.public import PublicKey, SealedBox
        pk = PublicKey(base64.b64decode(key_b64))
        sealed = SealedBox(pk).encrypt(plaintext.encode("utf-8"))
    except Exception:
        return None
    body = "\n".join(textwrap.wrap(base64.b64encode(sealed).decode("ascii"), 64))
    return f"{_ARMOR_HEADER}\n{body}\n{_ARMOR_FOOTER}"


def decrypt_blob(armored: str, private_key_b64: str) -> str:
    """Decrypt an armored block produced by :func:`encrypt_for_maintainer`.

    Tolerates surrounding text (e.g. the whole issue body pasted in): it reads
    only the content between the armor markers, and falls back to treating the
    entire input as the base64 body if no markers are present.
    """
    from nacl.public import PrivateKey, SealedBox

    body_lines, in_body, saw_markers = [], False, False
    for line in armored.splitlines():
        s = line.strip()
        if s == _ARMOR_HEADER:
            in_body, saw_markers = True, True
            continue
        if s == _ARMOR_FOOTER:
            break
        if in_body and s:
            body_lines.append(s)
    if not saw_markers:
        body_lines = [l.strip() for l in armored.splitlines()
                      if l.strip() and "```" not in l]
    if not body_lines:
        raise ValueError("No encrypted content found in the provided input.")

    sealed = base64.b64decode("".join(body_lines))
    sk = PrivateKey(base64.b64decode(private_key_b64.strip()))
    return SealedBox(sk).decrypt(sealed).decode("utf-8")


# ── Maintainer CLI entry points (called from main.py) ──────────────────────

def cli_generate_keypair() -> int:
    """Generate and store a feedback keypair; print the public key to embed."""
    from rich.console import Console
    console = Console()

    if not crypto_available():
        console.print("[red]PyNaCl is not installed.[/red] Install it with: "
                      "pip install pynacl", highlight=False)
        return 1

    if PRIVATE_KEY_PATH.exists():
        console.print(
            f"[yellow]A feedback private key already exists at "
            f"{PRIVATE_KEY_PATH}.[/yellow]", highlight=False)
        console.print(
            "[yellow]Refusing to overwrite it — that key decrypts existing "
            "reports. Remove it by hand first if you truly want a new one."
            "[/yellow]")
        return 1

    pub_b64, priv_b64 = generate_keypair()

    PRIVATE_KEY_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(PRIVATE_KEY_PATH, "w") as f:
        f.write(priv_b64 + "\n")
    try:
        os.chmod(PRIVATE_KEY_PATH, 0o600)
    except OSError:
        pass

    console.print("\n[bold green]✓ Feedback keypair generated.[/bold green]\n")
    console.print(f"Private key saved to: [bright_blue]{PRIVATE_KEY_PATH}[/bright_blue] "
                  f"(mode 600).", highlight=False)
    console.print("[yellow]Back this file up. If you lose it, encrypted reports "
                  "can never be read.[/yellow]\n")
    console.print("Embed this public key in "
                  "[bright_blue]src/proprep/utils/feedback_crypto.py[/bright_blue]:",
                  highlight=False)
    console.print(f'\n    FEEDBACK_PUBLIC_KEY_B64 = "{pub_b64}"\n', highlight=False)
    return 0


def cli_decrypt(source: str) -> int:
    """Decrypt an armored block from a file (or '-' for stdin) and print it."""
    import sys
    from rich.console import Console
    console = Console(stderr=True)

    if not crypto_available():
        console.print("[red]PyNaCl is not installed.[/red] Install it with: "
                      "pip install pynacl", highlight=False)
        return 1
    if not PRIVATE_KEY_PATH.exists():
        console.print(f"[red]No private key at {PRIVATE_KEY_PATH}.[/red] "
                      "Run: proprep --generate-feedback-keypair", highlight=False)
        return 1

    if source == "-":
        console.print("[grey50]Paste the encrypted block, then Ctrl-D:[/grey50]")
        armored = sys.stdin.read()
    else:
        try:
            armored = Path(source).read_text()
        except OSError as e:
            console.print(f"[red]Could not read {source}: {e}[/red]", highlight=False)
            return 1

    private_key_b64 = PRIVATE_KEY_PATH.read_text().strip()
    try:
        plaintext = decrypt_blob(armored, private_key_b64)
    except Exception as e:
        console.print(f"[red]Decryption failed: {e}[/red]", highlight=False)
        return 1

    # Plaintext transcript to stdout so it can be redirected to a file.
    print(plaintext)
    return 0
