import pyarrow as pa


def do_serialize_to_disk(
    data, filename, overwrite=True, skip_nones=True, format="pickle"
):
    if skip_nones:
        data = {k: v for k, v in data.items() if v is not None}
    if not overwrite and os.path.exists(filename):
        raise FileExistsError
    with open(filename, "wb") as f:
        if format == "arrow":
            buf = pa.serialize(data).to_buffer()
            f.write(buf)
        elif format == "pickle":
            pickle.dump(data, f)
    return data


def do_save_trenches(trenches, filename, overwrite=True):
    trenches = pd.concat(trenches)
    processing.write_dataframe_to_parquet(
        filename, trenches, merge=False, overwrite=overwrite
    )
    return trenches


def do_measure_and_write(
    trenches,
    frames,
    return_none=True,
    write=True,
    # filename_func=filename_func,
    **kwargs,
):
    if trenches is None:
        return None
    trenches = filter_trenches(trenches)
    res = measure(trenches, frames, **kwargs)
    if write:
        processing.write_images_and_measurements(
            res,
            filename_func=filename_func,
            dataframe_format="parquet",
            write_images=True,
            write_measurements=True,
        )
    if return_none:
        return None
    else:
        return res
