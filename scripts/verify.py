"""Verify the distribution and optionally replay its four finite checks."""

import argparse
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
COMPANION = ROOT / "companion"


def check_manifest():
    manifest = ROOT / "MANIFEST_SHA256.txt"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError("The distribution manifest is empty.")
    for line in lines:
        expected, relative_path = line.split("  ", 1)
        path = ROOT / relative_path
        if not path.resolve().is_relative_to(ROOT):
            raise ValueError(f"Path outside the distribution: {relative_path}")
        actual = sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch: {relative_path}")
    print(f"PASS manifest ({len(lines)} files)", flush=True)


def replay_examples():
    with TemporaryDirectory(prefix="replay-", dir=COMPANION) as temporary:
        output_dir = Path(temporary).relative_to(COMPANION)
        for example in ("m", "f"):
            reference = COMPANION / f"evidence/{example}-normal.json"
            for mode, flags in (("normal", []), ("optimized", ["-O"])):
                output = output_dir / f"{example}-{mode}.json"
                subprocess.run(
                    [
                        sys.executable,
                        *flags,
                        "-B",
                        f"scripts/check_example_{example}.py",
                        "--output",
                        str(output),
                    ],
                    cwd=COMPANION,
                    check=True,
                )
                if (COMPANION / output).read_bytes() != reference.read_bytes():
                    raise ValueError(f"Replay differs from reference: {example.upper()} {mode}")
                print(f"PASS {example.upper()} {mode} (reference bytes match)", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", action="store_true", help="run M and F in both Python modes")
    args = parser.parse_args()
    try:
        check_manifest()
        if args.replay:
            replay_examples()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
