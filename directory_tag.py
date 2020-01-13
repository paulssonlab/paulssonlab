#!/usr/bin/env python

import sys
import pathlib
import subprocess
from collections import Counter

MAX_TAGS = 3

ret = subprocess.run(
    ["git", "diff-index", "-z", "--cached", "HEAD", "--name-only"], capture_output=True
)

if ret.returncode != 0:
    sys.exit(ret.returncode)

files = [pathlib.Path(f.decode()) for f in ret.stdout.split(b"\x00") if f]
counts = Counter([f.parts[0] for f in files if len(f.parts) >= 2])

if not counts or len(counts) > MAX_TAGS:
    sys.exit(0)

tags = [e[0] for e in counts.most_common()]

tag_str = f"[{','.join(tags)}]"

print(tag_str)

with open(sys.argv[1], "r+") as f:
    msg = f.read()
    f.seek(0)
    f.write(tag_str + " " + msg)
