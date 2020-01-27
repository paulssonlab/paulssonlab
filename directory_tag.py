#!/usr/bin/env python

import sys
import pathlib
import subprocess
from collections import Counter

MAX_TAGS = 2

ret = subprocess.run(
    ["git", "diff-index", "-z", "--cached", "HEAD", "--name-only"], capture_output=True
)

if ret.returncode != 0:
    sys.exit(ret.returncode)

files = [pathlib.Path(f.decode()) for f in ret.stdout.split(b"\x00") if f]
counts = Counter([f.parent.parts[:2] for f in files])
del counts[()]
if len(counts) > MAX_TAGS:
    counts = Counter([f.parent.parts[:1] for f in files])
    del counts[()]
if not counts or len(counts) > MAX_TAGS:
    sys.exit(0)
tags = ["/".join(e) for e, _ in counts.most_common()]

tag_str = f"[{','.join(tags)}]"

with open(sys.argv[1], "r+") as f:
    msg = f.read()
    f.seek(0)
    f.write(tag_str + " " + msg)
