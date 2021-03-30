import pandas as pd
import re
from ast import literal_eval
from datetime import datetime
from cytoolz import get_in

POSITIONS_KEY = (b"SLxExperiment", b"ppNextLevelEx", b"", b"uLoopPars", b"Points", b"")


def get_grid(nd2):
    positions = get_in(POSITIONS_KEY, nd2._parser._raw_metadata.image_metadata)
    grid = [
        re.match(r"([^0-9]+)(\d+)", pos[b"dPosName"].decode()).groups()
        for pos in positions
    ]
    grid_df = pd.DataFrame(grid, columns=["row", "column"])
    grid_df.index.name = "pos"
    return grid_df


def parse_mux_log(filename):
    log = []
    with open(filename) as f:
        for line in f:
            m = re.match(r"^(\S+) on (\w+) sent (.*)$", line)
            time = datetime.fromisoformat(m.group(1))
            port = m.group(2)
            command = literal_eval(m.group(3))
            log.append((time, port, command))
    return log
