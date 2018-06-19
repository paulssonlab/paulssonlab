import numpy as np
import pandas as pd
from cytoolz import get_in


def get_position_metadata(metadata, grid_coords=True, reverse_grid="x"):
    def position_dataframe(d):
        df = pd.DataFrame.from_dict(d)
        df.rename(
            columns={
                "dPosName": "position_name",
                "dPosX": "x",
                "dPosY": "y",
                "dPosZ": "z",
                "dPFSOffset": "pfs_offset",
            },
            inplace=True,
        )
        df = df[["position_name", "x", "y", "z", "pfs_offset"]]
        if grid_coords:
            for coord in ("x", "y", "z"):
                coords = df[coord].unique()
                coords.sort()
                if coord in reverse_grid:
                    coords = coords[::-1]
                df[coord + "_idx"] = df[coord].map(
                    lambda c: np.where(coords == c)[0][0]
                )
        return df

    positions = pd.concat(
        {
            filename: position_dataframe(
                [
                    p
                    for p in get_in(
                        [
                            "image_metadata",
                            "SLxExperiment",
                            "ppNextLevelEx",
                            "",
                            "uLoopPars",
                            "Points",
                            "",
                        ],
                        md,
                    )
                ]
            )
            for filename, md in metadata.items()
        }
    )
    positions.index.names = ["filename", "position"]
    return positions
