import itertools as it
from collections.abc import Callable, Hashable, MutableMapping, MutableSequence
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any

import pandas as pd
from deltalake import write_deltalake

from paulssonlab.util.core import (
    ItemProxy,
    first,
    flatten_dict,
    iter_recursive,
    map_recursive,
)


class DelayedNotReadyError(RuntimeError):
    pass


class Delayed:
    pass


# TODO: use the sentinel library to give this a better __str__/__repr__
NotReady = object()


@dataclass(kw_only=True)
class DelayedValue(Delayed):
    value: Any = NotReady

    def is_ready(self):
        return self.value is not NotReady

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
    _result: Any = NotReady

    @cached_property
    def dependencies(self) -> "list[Delayed]":
        return list(get_delayed(self.func, self.args, self.kwargs))

    def is_ready(self):
        # print(f"IS READY ({len(self.dependencies)}):", all(dep.is_ready() for dep in self.dependencies))
        return self._result is not NotReady or all(
            dep.is_ready() for dep in self.dependencies
        )

    @property
    def is_pending(self):
        return self._result is NotReady

    def result(self):
        if (res := self._result) is NotReady:
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
        return res


@dataclass(kw_only=True)
class DelayedGetitem(Delayed):
    obj: MutableMapping | MutableSequence | Delayed
    key: Hashable | int | Delayed

    # @cached_property
    # def dependencies(self):
    #     return list(get_delayed(self.obj, self.key, recursive=False))

    def is_ready(self):
        if (isinstance(self.obj, Delayed) and not self.obj.is_ready()) or (
            isinstance(self.key, Delayed) and not self.key.is_ready()
        ):
            return False
        if self.key in self.obj:
            if isinstance(value := self.obj[self.key], Delayed):
                return value.is_ready()
            else:
                return True
        else:
            return False

    # @property
    # def value(self):
    #     return self.obj[self.key]

    # TODO: do we ever want to use this or will this only lead to confusion?
    # @value.setter
    # def value(self, new_value):
    #     self.obj[self.key] = new_value

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


def is_ready(value):
    if isinstance(value, Delayed):
        return value.is_ready()
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
        # print("POLLING")
        ids = list(self._items.keys())
        while True:
            any_fired = False
            idx = 0
            while idx < len(ids):
                # print(f"IDX: {idx}/{len(ids)}")
                id_ = ids[idx]
                delayed = self._items[id_]
                # print("*", delayed, "IS_READY", delayed.is_ready())
                if delayed.is_ready():
                    # print("!!! FIRING")
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
                        # print(f"ADDING DEP {delayed} -> {dep}")
                        self.append(dep)
            self._items[id(delayed)] = delayed


class DelayedStore:
    def __init__(self, queue):
        self.value = {}
        self.queue = queue
        self._write_queue = set()
        self.delayed = ItemProxy(self, "delayed")

    def __len__(self):
        return len(self.value)

    def __getitem__(self, key):
        if key in self:
            return get_result_if_ready(self.value[key])
        else:
            return DelayedGetitem(obj=self, key=key)

    def _delayed_getitem(self, key):
        if key in self:
            return self.value[key]
        else:
            return DelayedGetitem(obj=self, key=key)

    def __setitem__(self, key, value):
        if key in self.value:
            raise RuntimeError(f"store is immutable, cannot overwrite key {key}")
        if isinstance(value, DelayedCallable):
            self.queue.append(value)
        self._write_queue.add(key)
        self.value[key] = value

    def _delayed_setitem(self, key, value):
        self.__setitem__(key, value)

    def __delitem__(self, key):
        raise RuntimeError(f"store is immutable, cannot delete key {key}")

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
        return self.value.clear()

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
        return {k: get_result_if_ready(v) for k, v in self.value.items()}.items()

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
        if key in self.value:
            # print("A", key, value)
            return self[key]
        else:
            # print("B", key, value)
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


class DelayedArrayStore(DelayedStore):
    def __init__(self, queue, output_path):
        self.output_path = output_path
        super().__init__(queue)


class DelayedTableStore(DelayedStore):
    DEFAULT_WRITE_OPTIONS = {"mode": "append", "engine": "rust"}

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
        self.queue.append(writer)
        self._write_queue.clear()

    @staticmethod
    def _concat_table(items, schema):
        # TODO: replace pd.concat with faster arrow operation?
        items = flatten_dict(items, concatenate_keys=True)
        first_key = first(items)
        names = [s[0] for s in schema[: len(first_key)]]
        return pd.concat(items, names=names)

    @classmethod
    def _write(cls, items, schema, path, write_options):
        table = cls._concat_table(items, schema)
        write_deltalake(path, table, **write_options)
