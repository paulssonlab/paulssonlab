import numpy
import pandas
import tables

"""
Utility functions for taking a Pandas dataframe & writing it to HDF5.
For whatever reason, the to_hdf() functionality doesn't seem to write to the exact same format
as PyTables tables written natively, so I've re-implented the feature.
This seems to be because Pandas either stores strings as generic objects,
or as a custom StringDtype.

The solution is to convert text columns to utf-8 when working with them,
and then back to numpy.bytes_ when writing them to disk.
Actually, I'm a bit lost as to when it's necessary to convert to utf-8,
but numpy.bytes_ are for sure required at the very end.

And, we have to define the table's list of columns, and dynamically
convert between dataframe columns (with numpy datatypes) and pytables column types.
"""


def dtypes_to_definitions(name, dtypes):
    """
    Input a dataframe dtype specification as a dict,
    and a name for the type of data,
    return a pytables row definition as a dict.
    """
    # FIXME Other missing dtypes? To be completed as needed.
    # Special case for strings - see below.
    conversion_dict = {
        numpy.dtype("bool"): tables.BoolCol(),
        numpy.dtype("int8"): tables.Int8Col(),
        numpy.dtype("int16"): tables.Int16Col(),
        numpy.dtype("int32"): tables.Int32Col(),
        numpy.dtype("int64"): tables.Int64Col(),
        numpy.dtype("uint8"): tables.UInt8Col(),
        numpy.dtype("uint16"): tables.UInt16Col(),
        numpy.dtype("uint32"): tables.UInt32Col(),
        numpy.dtype("uint64"): tables.UInt64Col(),
        numpy.dtype("float32"): tables.Float32Col(),
        numpy.dtype("float64"): tables.Float64Col(),
    }

    new_dtypes = {}

    for label, dtype in dtypes.items():
        # Is it a string stored as numpy bytes array?
        if dtype.kind == "S":
            new_dtypes[label] = tables.StringCol(dtype.itemsize)
        else:
            new_dtypes[label] = conversion_dict[dtype]

    return type(name, (tables.IsDescription,), new_dtypes)


def write_dataframe_to_hdf5(dataframe, outfile, name):
    """
    Given a pandas dataframe, and an out path, and a name for the file,
    write the dataframe to HDF5 using a more standard pytables format
    than provided by pandas' to_hdf() function.
    """
    dtypes_dict = dataframe.dtypes.to_dict()
    row_definition = dtypes_to_definitions(name, dtypes_dict)
    h5file = tables.open_file(outfile, "w")

    # Don't bother including a title (optional argument)
    t = h5file.create_table(
        "/", name, row_definition, filters=tables.Filters(complevel=1, complib="zlib")
    )
    data_row = t.row

    # All the rows
    for entry in dataframe.itertuples():
        # All the columns
        for column in dtypes_dict.keys():
            # getattr allows us to get the attribute by variable
            data_row[column] = getattr(entry, column)
        # Write the data
        data_row.append()

    # Done
    t.flush()
    t.close()
    h5file.close()
