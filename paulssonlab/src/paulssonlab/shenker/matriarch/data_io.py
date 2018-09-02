import pyarrow as pa
import pyarrow.parquet as pq
import numpy_indexed as npi
from cytoolz import take
import json
from util import grouper, tqdm_auto

DEFAULT_ROW_GROUP_SIZE = 1_000_000


def read_parquet(
    source,
    columns=None,
    categories=None,
    length=None,
    strings_to_categories=True,
    progress_bar=tqdm_auto,
):
    reader = pq.ParquetFile(source)
    tables = []
    category_idxs = None
    if progress_bar is not None:
        pbar = progress_bar(total=reader.metadata.num_rows)
    row_groups = range(length or reader.num_row_groups)
    for row_group in row_groups:
        table = reader.read_row_group(row_group, nthreads=4, columns=columns)
        if progress_bar is not None:
            pbar.update(len(table))
        if strings_to_categories:
            new_categories = []
            for name in table.schema.names:
                if table.schema.field_by_name(name).type == pa.string():
                    new_categories.append(name)
            categories = list(set((categories or []) + new_categories))
        if categories is not None:
            if category_idxs is None:
                category_idxs = [
                    table.schema.get_field_index(column) for column in categories
                ]
            for idx in category_idxs:
                table = table.set_column(idx, table.column(idx).dictionary_encode())
        tables.append(table)
    if categories is not None:
        for idx in category_idxs:
            category_columns = [table.column(idx) for table in tables]
            new_category_columns = harmonize_dictionaries(category_columns)
            for i, new_column in enumerate(new_category_columns):
                tables[i] = tables[i].set_column(idx, new_column)
    if progress_bar is not None:
        pbar.close()
    return pa.concat_tables(tables)


def first_index(batch):
    index_columns = json.loads(batch.schema.metadata[b"pandas"])["index_columns"]
    return tuple(
        batch.column(batch.schema.get_field_index(name))[0].as_py()
        for name in index_columns
    )


def sort_arrow_to_parquet(
    arrow_filename,
    parquet_filename=None,
    length=None,
    row_group_size=DEFAULT_ROW_GROUP_SIZE,
    progress_bar=tqdm_auto,
    sort_within_batch=False,
):
    if parquet_filename is None:
        parquet_filename = arrow_filename.replace(".arrow", ".parquet")
    arrow_file = pa.OSFile(arrow_filename)
    import io

    arrow_file = io.BufferedReader(arrow_file, io.DEFAULT_BUFFER_SIZE * 10000)
    # arrow_file = pa.memory_map(arrow_filename)
    reader = pa.open_stream(arrow_file)
    if progress_bar is not None:
        pbar = progress_bar(desc="reading")
    if length is not None:
        reader = take(length, reader)
    batches = []
    for batch in reader:
        if progress_bar is not None:
            pbar.update(len(batch))
        if sort_within_batch:
            raise NotImplementedError
        batches.append((first_index(batch), batch))
    table = pa.Table.from_batches([b[1] for b in sorted(batches)])
    with pq.ParquetWriter(parquet_filename, table.schema) as writer:
        writer.write_table(table, row_group_size=row_group_size)
    arrow_file.close()


def copy_arrow2(
    in_filename, out_filename, length=None, batch_size=10000, process_func=None
):
    in_file = pa.memory_map(in_filename)
    reader = pa.RecordBatchStreamReader(in_file)
    out_file = pa.OSFile(out_filename, "wb")
    table0 = pa.Table.from_batches([reader.read_next_batch()])
    if process_func is not None:
        table0 = process_func(table0)
    writer = pa.RecordBatchStreamWriter(out_file, table0.schema)
    writer.write_table(table0)
    if length is not None:
        reader = take(length, reader)
    t0 = time.time()
    for i, batches in enumerate(grouper(reader, batch_size)):
        if True:  # i % 100 == 0:
            t = time.time()
            dt = t - t0
            t0 = t
            print("batch", i, "time {:.2f}".format(dt))
        table = pa.Table.from_batches(batches)  # .drop(columns_to_drop)
        if process_func is not None:
            table = process_func(table)
        print("    rows per second", len(table) / dt)
        writer.write_table(table)


def read_arrow(arrow_filename, columns, categorical_columns=None, batch_size=1000):
    reader = pa.open_stream(arrow_filename)
    table0 = pa.Table.from_batches([reader.read_next_batch()])
    columns_to_drop = list(set(c.name for c in table0.columns) - set(columns))
    table1 = table0.drop(columns_to_drop)
    imos = pa.BufferOutputStream()
    writer = pa.RecordBatchStreamWriter(imos, table1.schema)
    t0 = time.time()
    for i, batches in enumerate(grouper(reader, batch_size)):
        if True:  # i % 100 == 0:
            t = time.time()
            dt = t - t0
            t0 = t
            print("batch", i, "time {:.2f}".format(dt))
        table = pa.Table.from_batches(batches).drop(columns_to_drop)
        for i in range(table.num_columns):
            if table.column(i).name == "filename":
                table = table.set_column(i, table.column(i).dictionary_encode())
        print("    rows per second", len(table) / dt)
        writer.write_table(table)
    # with pq.ParquetWriter(parquet_filename, table0.schema) as writer:
    #     writer.write_table(table0)
    #    for batches in util.grouper(reader, batch_size):
    #         table = pa.Table.from_batches(batches)
    #         writer.write_table(table)
    output_reader = pa.open_stream(imos.getvalue())
    return output_reader


def harmonize_dictionaries(columns):
    dictionaries = [
        list(map(str, column.data.chunk(0).dictionary)) for column in columns
    ]
    # sorting is important
    # otherwise sorting by string will not lexsort pandas categoricals
    master_dictionary = pa.array(sorted(set.union(*(set(d) for d in dictionaries))))
    value_to_index = dict(map(reversed, enumerate(master_dictionary)))
    new_columns = []
    for column in columns:
        column_name = column.name
        rewrite_rules = {
            orig: value_to_index[value]
            for orig, value in enumerate(column.data.chunk(0).dictionary)
        }
        from_values = list(rewrite_rules.keys())
        to_values = list(rewrite_rules.values())
        new_chunks = []
        for chunk_idx in range(column.data.num_chunks):
            chunk = column.data.chunk(chunk_idx)
            ary = chunk.indices.to_numpy()
            npi.remap(ary, from_values, to_values, inplace=True)
            new_chunk = pa.DictionaryArray.from_arrays(chunk.indices, master_dictionary)
            new_chunks.append(new_chunk)
        chunked_ary = pa.chunked_array(new_chunks)
        new_column = pa.Column.from_array(column_name, chunked_ary)
        new_columns.append(new_column)
    return new_columns


def arrow_to_parquet(
    arrow_filename, parquet_filename=None, batch_size=1000, length=None
):
    if parquet_filename is None:
        parquet_filename = arrow_filename.replace(".arrow", ".parquet")
    arrow_file = pa.OSFile(arrow_filename)
    # arrow_file = pa.memory_map(arrow_filename)
    # parquet_mmap = pa.memory_map(parquet_filename, 'wb')
    reader = pa.open_stream(arrow_file)
    table0 = pa.Table.from_batches([reader.read_next_batch()])
    if length is not None:
        reader = take(length, reader)
    with pq.ParquetWriter(parquet_filename, table0.schema) as writer:
        # with pq.ParquetWriter(parquet_mmap, table0.schema) as writer:
        writer.write_table(table0)
        for batches in grouper(reader, batch_size):
            table = pa.Table.from_batches(batches)
            writer.write_table(table)
    # arrow_file.close()
    # parquet_mmap.flush()
    # parquet_mmap.close()
