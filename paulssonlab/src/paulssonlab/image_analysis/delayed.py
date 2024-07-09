import itertools as it
from collections import defaultdict
from collections.abc import Callable, Hashable, MutableMapping, MutableSequence
from dataclasses import dataclass
from functools import cached_property
from numbers import Integral
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr
from cytoolz import merge

from paulssonlab.image_analysis.image import pad
from paulssonlab.util.core import (
    ItemProxy,
    extract_keys,
    first,
    flatten_dict,
    iter_recursive,
    map_recursive,
    str_placeholders,
)


class DelayedNotReadyError(RuntimeError):
    pass


class Delayed:
    pass


# TODO: use the sentinel library to give this a better __str__/__repr__
NotAvailable = object()


@dataclass(kw_only=True)
class DelayedValue(Delayed):
    value: Any = NotAvailable

    def is_ready(self):
        return self.value is not NotAvailable

    is_available = is_ready

    def result(self):
        if not self.is_ready():
            raise DelayedNotReadyError
        return self.value


@dataclass(kw_only=True)
class DelayedCallable(Delayed):
    func: Callable | None
    args: list | None = None
    kwargs: dict | None = None
    _keep_args: bool = False
    _keep_dependencies: bool = False
    _result: Any = NotAvailable

    @cached_property
    def dependencies(self) -> "list[Delayed]":
        return list(get_delayed(self.func, self.args, self.kwargs))

    def is_ready(self):
        return self.is_available() or all(dep.is_ready() for dep in self.dependencies)

    def is_available(self):
        return self._result is not NotAvailable

    def result(self):
        if (res := self._result) is NotAvailable:
            func = unbox_delayed(self.func)
            args = unbox_delayed(self.args) or ()
            kwargs = unbox_delayed(self.kwargs) or {}
            if self._keep_args:
                # keep unboxed values for debugging
                # (e.g., to make it easy to re-run function call)
                # that's why we do this before running function call, in case that fails
                self.func = func
                self.args = args
                self.kwargs = kwargs
            res = func(*args, **kwargs)
            self._result = res
            if not self._keep_args:
                self.func = None
                self.args = None
                self.kwargs = None
            if not self._keep_dependencies:
                # keep dependencies for debugging
                # note that this will prevent garbage-collecting unused results of dependencies
                # (e.g., dask Futures)
                self.dependencies = []
        return res


@dataclass(kw_only=True)
class DelayedGetitem(Delayed):
    obj: MutableMapping | MutableSequence
    key: Hashable | int

    def is_ready(self):
        if self.key in self.obj:
            if isinstance(value := self.obj[self.key], Delayed):
                return value.is_ready()
            else:
                return True
        else:
            return False

    def is_available(self):
        if self.key in self.obj:
            if isinstance(value := self.obj[self.key], Delayed):
                return value.is_available()
            else:
                return True
        else:
            return False

    def result(self):
        if not self.is_ready():
            raise DelayedNotReadyError
        obj = get_result(self.obj)
        key = get_result(self.key)
        value = obj[key]
        return get_result(value)


def get_delayed(*args, recursive=True):
    # TODO: look how dask implements unboxing
    if recursive:
        args = it.chain.from_iterable((iter_recursive(arg) for arg in args))
    return filter(lambda x: isinstance(x, Delayed), args)


def get_result(value):
    if isinstance(value, Delayed):
        return value.result()
    else:
        return value


def get_result_if_ready(value):
    if isinstance(value, Delayed) and value.is_ready():
        return value.result()
    else:
        return value


def get_result_if_available(value):
    if isinstance(value, Delayed) and value.is_available():
        return value.result()
    else:
        return value


def is_ready(value):
    if isinstance(value, Delayed):
        return value.is_ready()
    else:
        return True


def is_available(value):
    if isinstance(value, Delayed):
        return value.is_available()
    else:
        return True


def unbox_delayed(obj):
    return map_recursive(get_result, obj)


class DelayedQueue:
    def __init__(self, wrapper=None):
        self._items = {}
        self.wrapper = wrapper

    def delayed(self, func, *args, **kwargs):
        _keep_args = kwargs.pop("_keep_args", False)
        wrapper = kwargs.pop("_wrapper", self.wrapper)
        if wrapper is not None:
            func = self.wrapper(func)
        dc = DelayedCallable(
            func=func,
            args=args,
            kwargs=kwargs,
            _keep_args=_keep_args,
        )
        return dc

    def poll(self):
        ids = list(self._items.keys())
        while True:
            any_fired = False
            idx = 0
            while idx < len(ids):
                id_ = ids[idx]
                delayed = self._items[id_]
                if delayed.is_ready():
                    delayed.result()
                    self._items.pop(id_, None)
                    del ids[idx]
                    any_fired = True
                else:
                    idx += 1

            if not any_fired:
                break
        return

    def append(self, delayed):
        if id(delayed) not in self._items:
            if hasattr(delayed, "dependencies") and delayed.dependencies:
                for dep in delayed.dependencies:
                    if isinstance(dep, DelayedCallable):
                        self.append(dep)
            self._items[id(delayed)] = delayed


class DelayedStore:
    def __init__(self, queue):
        self.value = {}
        self.queue = queue
        self._write_queue = set()
        self.writers = {}
        self._evicted = set()
        self.delayed = ItemProxy(self, "delayed")

    def __len__(self):
        return len(self.value)

    def __getitem__(self, key):
        if key in self:
            return get_result_if_available(self.value[key])
        else:
            return DelayedGetitem(obj=self, key=key)

    def _delayed_getitem(self, key):
        if key in self:
            return self.value[key]
        else:
            return DelayedGetitem(obj=self, key=key)

    def __setitem__(self, key, value):
        if key in self.value or key in self._evicted:
            raise RuntimeError(f"store is immutable, cannot overwrite key {key}")
        if isinstance(value, DelayedCallable):
            self.queue.append(value)
        self._write_queue.add(key)
        self.value[key] = value

    def _delayed_setitem(self, key, value):
        self.__setitem__(key, value)

    def __delitem__(self, key):
        self._evicted.add(key)
        if key not in self._write_queue:
            self.value.pop(key, None)

    def _delayed_delitem(self, key):
        return self.__delitem__(self, key)

    def __contains__(self, key):
        return key in self.value

    def _delayed_contains(self, key):
        return self.__contains__(key)

    def __iter__(self):
        return iter(self.value)

    def _delayed_iter(self):
        return self.__iter__()

    def clear(self):
        self.value.clear()
        self._write_queue.clear()

    def __copy__(self):
        new = type(self).__new__()
        new.value = self.value.copy()
        return new

    def copy(self):
        return self.__copy__()

    @classmethod
    def fromkeys(cls, iterable, value):
        raise NotImplementedError

    def get(self, key, default):
        if key in self:
            return self[key]
        else:
            return default

    def items(self):
        # TODO: this ItemsView does not reference underlying dict self.value
        # (i.e., has correct type but doesn't update when new keys are added)
        # this shouldn't be a big deal though
        return {k: get_result_if_available(v) for k, v in self.value.items()}.items()

    def _delayed_items(self):
        return self.value.items()

    def keys(self):
        return self.value.keys()

    def pop(self, key, default=None):
        raise NotImplementedError

    def popitem(self):
        raise NotImplementedError

    def __reversed__(self):
        return reversed(self.value)

    def setdefault(self, key, value):
        if key in self:
            return get_result_if_available(self.value[key])
        else:
            self[key] = value
            return value

    def update(self, *args, **kwargs):
        raise NotImplementedError

    def values(self):
        return self.value.values()

    def __or__(self, other):
        raise NotImplementedError

    def __ior__(self, other):
        raise NotImplementedError

    def write(self):
        raise NotImplementedError

    def flush(self):
        return self.write()


def _shape_for_subarrays(shapes, pad_if_needed=True):
    if pad_if_needed:
        shapes = np.asarray(shapes)
        return tuple(int(x) for x in np.max(shapes, axis=0))
    else:
        if not all(shape == shapes[0] for shape in shapes[1:]):
            raise ValueError(f"got heterogenous shapes: {list(set(shapes))}")
        return shapes[0]


def _dtype_for_arrays(dtypes, fill_value=None):
    if not all(dtype == dtypes[0] for dtype in dtypes[1:]):
        raise ValueError(f"got heterogenous dtypes: {list(set(dtypes))}")
    dtype = dtypes[0]
    if fill_value is not None:
        dtype = np.promote_types(dtype, np.min_scalar_type(fill_value))
    return dtype


def _shape_for_keys(keys, old_shape=None):
    keys = np.asarray(keys)
    if old_shape is not None:
        keys = np.concatenate((keys, np.asarray(old_shape)[np.newaxis, :] - 1))
    return tuple(int(x) for x in (np.max(keys, axis=0) + 1))


class DelayedZarrStore(DelayedStore):
    DEFAULT_WRITE_OPTIONS = {"fill_value": np.nan}

    def __init__(self, queue, output_path, write_options=None):
        self.output_path = output_path
        self.write_options = {**self.DEFAULT_WRITE_OPTIONS, **(write_options or {})}
        super().__init__(queue)

    def write(self):
        # we directly access self.value[k] instead of self[k] so we don't check readiness
        items_to_write = {
            k: self.value[k] for k in self._write_queue if k in self.value
        }
        if not items_to_write:
            return
        writer = self.queue.delayed(
            self._write,
            items_to_write,
            self.output_path,
            self.write_options,
        )
        self.writers[tuple(items_to_write.keys())] = writer
        self.queue.append(writer)
        self._write_queue -= items_to_write.keys()

    @staticmethod
    def _write(items, path, write_options):
        chunks_spec = write_options.get("chunks", None)
        fill_value = write_options.get("fill_value", None)
        if chunks_spec is None:
            chunks_spec = ()
        for store_key, arrays in items.items():
            if (num_placeholders := len(str_placeholders(str(path)))) != len(store_key):
                raise ValueError(
                    f"expected {len(store_key)} placeholders in output_path, got {num_placeholders} instead"
                )
            path = Path(str(path).format(*store_key))
            store = zarr.DirectoryStore(path)
            arrays = flatten_dict(arrays, concatenate_keys=True)
            array_shape = tuple(
                np.max([array.shape for array in arrays.values()], axis=0)
            )
            if not path.exists():
                key_shape = tuple(np.max(list(arrays.keys()), axis=0) + 1)
                dtype = _dtype_for_arrays(
                    [array.dtype for array in arrays.values()], fill_value=fill_value
                )
                shape = (*key_shape, *array_shape)
                chunks = chunks_spec + (None,) * (len(shape) - len(chunks_spec))
                ary = zarr.open_array(
                    store,
                    mode="a",
                    zarr_version=2,
                    shape=shape,
                    chunks=chunks,
                    dtype=dtype,
                    fill_value=fill_value,
                )
            else:
                ary = zarr.open_array(store, mode="a", zarr_version=2)
            for array_key, array in arrays.items():
                ary[array_key] = pad(array, array_shape, fill_value=fill_value)


# def _normalize_chunks(chunks, shape):
#     return tuple(c if c is not None else s for c, s in zip(chunks, shape))


def group_by_chunks(idxs, chunks):
    def chunk_index(idx):
        return tuple(i // c for i, c in zip(idx, chunks))

    groups = defaultdict(list)
    for key in idxs:
        groups[chunk_index(key)].append(key)
    return groups


def indices_for_chunk(chunk_idx, chunks):
    return set(
        it.product(*[range(i * c, (i + 1) * c) for i, c in zip(chunk_idx, chunks)])
    )


def complete_chunks(idxs, chunks):
    chunk_groups = group_by_chunks(idxs, chunks)
    return {
        k: v for k, v in chunk_groups.items() if set(v) == indices_for_chunk(k, chunks)
    }


def stack_subarrays(
    source_arrays,
    key_shape,
    subarray_shape,
    origin=None,
    fill_value=np.nan,
    dtype=None,
    pad_if_needed=True,
):
    if origin is None:
        origin = tuple(0 for _ in range(len(key_shape)))
    fill_array = np.full(subarray_shape, fill_value, dtype=dtype)
    sorted_keys = [
        key if key in source_arrays else None
        for key in it.product(*[range(o, o + s) for o, s in zip(origin, key_shape)])
    ]
    unused_keys = set(source_arrays.keys()) - set(sorted_keys)
    if unused_keys:
        raise ValueError(f"got extraneous keys: {unused_keys}")
    sorted_arrays = [
        fill_array if key is None else source_arrays[key].astype(dtype)
        for key in sorted_keys
    ]
    if pad_if_needed:
        sorted_arrays = [
            pad(array, subarray_shape, fill_value=fill_value) for array in sorted_arrays
        ]
    return np.stack(sorted_arrays).reshape(*key_shape, *subarray_shape)


class DelayedBatchedZarrStore(DelayedZarrStore):
    DEFAULT_WRITE_OPTIONS = {
        "pad_if_needed": True,
        "write_complete_chunks": True,
        "fill_value": np.nan,
        "store": zarr.DirectoryStore,
        "synchronizer": zarr.ProcessSynchronizer,
    }

    def __init__(self, queue, output_path, write_options=None):
        super().__init__(queue, output_path, write_options=write_options)
        chunks_spec = self.write_options["chunks"]
        if chunks_spec is None or any(not isinstance(x, Integral) for x in chunks_spec):
            raise ValueError(f"chunks must be a tuple of integers: {chunks_spec}")

    def write(self, partial_chunks=False):
        path_placeholders = str_placeholders(str(self.output_path))
        subarrays_by_path = defaultdict(dict)
        subarray_key_to_full_key = {}
        for key in self._write_queue:
            if key not in self.value:
                continue

            key_subset = tuple(
                key[idx] for idx in range(len(key)) if idx not in path_placeholders
            )
            remaining_placeholders = [p for p in path_placeholders if p >= len(key)]
            if remaining_placeholders:
                min_remaining_placeholder = min(remaining_placeholders)
                placeholder_replacements = [
                    f"{{{p-min_remaining_placeholder}}}" for p in remaining_placeholders
                ]
            else:
                placeholder_replacements = []
            path = str(self.output_path).format(*key, *placeholder_replacements)
            subarrays_by_path[path][key_subset] = self.value[key]
            subarray_key_to_full_key[(path, key_subset)] = key
        if partial_chunks:

            def chunk_filter_func(keys, chunks):
                return {None: keys}

        else:
            chunk_filter_func = complete_chunks
        items_to_write = {
            path: merge(
                [
                    extract_keys(subarrays, subarray_keys)
                    for _, subarray_keys in chunk_filter_func(
                        list(subarrays.keys()), self.write_options["chunks"]
                    ).items()
                ]
            )
            for path, subarrays in subarrays_by_path.items()
        }
        # filter out paths with no subarrays
        items_to_write = {
            path: subarrays for path, subarrays in items_to_write.items() if subarrays
        }
        if not items_to_write:
            return

        writer = self.queue.delayed(
            self._write,
            items_to_write,
            {"batch_chunks": not partial_chunks, **self.write_options},
        )
        # writer_key: ((path1, (chunk_key1.1, chunk_key1.2, ...)), (path2, (chunk_key2.1, ...)), ...)
        writer_key = tuple(
            (path, tuple(sorted(subarrays.keys())))
            for path, subarrays in sorted(items_to_write.items())
        )
        self.writers[writer_key] = writer
        self.queue.append(writer)
        keys_to_write = set(
            subarray_key_to_full_key[(path, subarray_key)]
            for path, subarrays in items_to_write.items()
            for subarray_key in subarrays.keys()
        )
        self._write_queue -= keys_to_write
        for key in keys_to_write & self._evicted:
            self.value.pop(key, None)

    def flush(self):
        return self.write(partial_chunks=True)

    @staticmethod
    def _write(items, write_options):
        pad_if_needed = write_options.get("pad_if_needed", True)
        batch_chunks = write_options.get("batch_chunks", False)
        chunks_spec = write_options.get("chunks", None)
        fill_value = write_options.get("fill_value", None)
        if chunks_spec is None:
            chunks_spec = ()
        items_by_path = defaultdict(dict)
        for path_template, path_subarrays in items.items():
            path_placeholders = str_placeholders(path_template)
            for key, value in path_subarrays.items():
                if isinstance(value, dict):
                    value = flatten_dict(value, concatenate_keys=True)
                else:
                    value = {(): value}
                for value_key, array in value.items():
                    if path_placeholders and (
                        (max_placeholder := max(path_placeholders)) + 1 > len(value_key)
                    ):
                        raise ValueError(
                            f"path '{path_template}' requires value key with length at least {max_placeholder+1}, but instead got value key: {value_key}"
                        )
                    value_key_subset = tuple(
                        value_key[idx]
                        for idx in range(len(value_key))
                        if idx not in path_placeholders
                    )
                    path = path_template.format(*value_key)
                    items_by_path[path][key + value_key_subset] = np.asarray(array)
        for store_path, subarrays in items_by_path.items():
            store = write_options["store"](store_path)
            synchronizer = write_options["synchronizer"](store_path + ".lock")
            all_subarrays = [subarray for subarray in subarrays.values()]
            all_keys = [key for key in subarrays.keys()]
            subarray_shape = _shape_for_subarrays(
                [subarray.shape for subarray in all_subarrays],
                pad_if_needed=pad_if_needed,
            )
            dtype = _dtype_for_arrays(
                [subarray.dtype for subarray in all_subarrays],
                fill_value=fill_value,
            )
            if not Path(store_path).exists():
                key_shape = _shape_for_keys(all_keys)
                new_shape = (*key_shape, *subarray_shape)
                chunks = chunks_spec + (None,) * (len(new_shape) - len(chunks_spec))
                dest_array = zarr.open_array(
                    store,
                    mode="a",
                    zarr_version=2,
                    shape=new_shape,
                    chunks=chunks,
                    dtype=dtype,
                    fill_value=fill_value,
                    synchronizer=synchronizer,
                )
            else:
                dest_array = zarr.open_array(
                    store, mode="a", zarr_version=2, synchronizer=synchronizer
                )
                key_shape = _shape_for_keys(
                    all_keys,
                    old_shape=dest_array.shape[
                        : len(dest_array.shape) - len(subarray_shape)
                    ],
                )
                new_shape = (
                    *key_shape,
                    *subarray_shape,
                )
            if dest_array.shape != new_shape:
                dest_array.resize(*new_shape)
            if batch_chunks:
                subarrays_by_chunk = {
                    chunk_key: extract_keys(subarrays, subarray_keys)
                    for chunk_key, subarray_keys in group_by_chunks(
                        list(subarrays.keys()),
                        chunks_spec,
                    ).items()
                }
                for chunk_key, chunk_subarrays in subarrays_by_chunk.items():
                    origin = tuple(k * c for k, c in zip(chunk_key, chunks_spec))
                    chunk_shape = tuple(
                        min(c, k - o) for c, k, o in zip(chunks_spec, key_shape, origin)
                    )
                    chunk_array = stack_subarrays(
                        chunk_subarrays,
                        chunk_shape,
                        subarray_shape,
                        origin=origin,
                        fill_value=fill_value,
                        dtype=dtype,
                        pad_if_needed=pad_if_needed,
                    )
                    dest_array.blocks[chunk_key] = chunk_array
            else:
                for key, subarray in subarrays.items():
                    print()
                    print(">>>>>>", key)
                    print("PRE-PAD", subarray.shape, subarray.dtype, subarray)
                    if pad_if_needed:
                        subarray = pad(
                            subarray.astype(dtype),
                            subarray_shape,
                            fill_value=fill_value,
                        )
                    print("POST-PAD", subarray.shape, subarray.dtype, subarray)
                    dest_array[key] = subarray
                    print("DONE SETTING")


class DelayedTableStore(DelayedStore):
    # DEFAULT_WRITE_OPTIONS = {"mode": "append", "engine": "rust"}
    DEFAULT_WRITE_OPTIONS = {"existing_data_behavior": "delete_matching"}

    def __init__(self, queue, output_path, schema=None, write_options=None):
        self.output_path = output_path
        self.schema = schema
        self.write_options = {**self.DEFAULT_WRITE_OPTIONS, **(write_options or {})}
        super().__init__(queue)

    def write(self):
        # we directly access self.value[k] instead of self[k] so we don't check readiness
        items_to_write = {
            k: self.value[k] for k in self._write_queue if k in self.value
        }
        if not items_to_write:
            return
        writer = self.queue.delayed(
            self._write,
            items_to_write,
            self.schema,
            self.output_path,
            self.write_options,
        )
        self.writers[tuple(items_to_write.keys())] = writer
        self.queue.append(writer)
        self._write_queue -= items_to_write.keys()

    @staticmethod
    def _concat_table(items, schema):
        # TODO: replace pd.concat with faster arrow operation?
        items = flatten_dict(items, concatenate_keys=True)
        first_key = first(items)
        names = [s[0] for s in schema[: len(first_key)]]
        df = pd.concat(items, names=names)
        if isinstance(df.index, pd.MultiIndex):
            levels_to_drop = [
                idx for idx, name in enumerate(df.index.names) if name is None
            ]
            if len(levels_to_drop) == len(df.index.names):
                levels_to_drop = levels_to_drop[1:]
            df = df.droplevel(level=levels_to_drop, axis=0)
            if df.index.names[-1] is not None:
                df = df.reset_index()
        else:
            if df.index.name is not None:
                df = df.reset_index()
        return df

    @classmethod
    def _write(cls, items, schema, path, write_options):
        df = cls._concat_table(items, schema)
        # write_deltalake(path, table, **write_options)
        pq.write_to_dataset(
            pa.Table.from_pandas(df, preserve_index=False), path, **write_options
        )
