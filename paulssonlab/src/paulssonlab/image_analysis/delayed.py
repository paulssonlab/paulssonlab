import itertools as it
from collections.abc import Callable, Hashable, MutableMapping, MutableSequence
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any

from paulssonlab.util.core import iter_recursive, map_recursive


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

    def result(self, **kwargs):
        if not self.is_ready():
            raise DelayedNotReadyError
        return self.value


@dataclass(kw_only=True)
class DelayedCallable(Delayed):
    func: Callable | None
    args: list | None = None
    kwargs: dict | None = None
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

    def result(self, clear_callable=True, **kwargs):
        if (res := self._result) is NotReady:
            func = unbox_delayed(self.func)
            args = unbox_delayed(self.args) or ()
            kwargs = unbox_delayed(self.kwargs) or {}
            print()
            print("&&&", func, args, kwargs)
            print()
            res = func(*args, **kwargs)
            self._result = res
            if clear_callable:
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

    def result(self, **kwargs):
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


def unbox_delayed(obj):
    return map_recursive(get_result, obj)


# def unbox_delayed(obj):
#     return map_recursive(
#         lambda x: x.result() if isinstance(x, Delayed) else x,
#         obj,
#     )


class DelayedQueue:
    def __init__(self):
        self._items = {}

    def delayed(self, func, *args, **kwargs):
        dc = DelayedCallable(
            func=func,
            args=args,
            kwargs=kwargs,
        )
        return dc

    def poll(self):
        print("POLLING")
        ids = list(self._items.keys())
        while True:
            any_fired = False
            idx = 0
            while idx < len(ids):
                print(f"IDX: {idx}/{len(ids)}")
                id_ = ids[idx]
                delayed = self._items[id_]
                print("*", delayed, "IS_READY", delayed.is_ready())
                if delayed.is_ready():
                    print("!!! FIRING")
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
                        print(f"ADDING DEP {delayed} -> {dep}")
                        self.append(dep)
            self._items[id(delayed)] = delayed


class DelayedStore:
    def __init__(self, queue):
        self.value = {}
        self.queue = queue

    def __len__(self):
        return len(self.value)

    def __getitem__(self, key):
        if key in self:
            value = self.value[key]
            return value
        else:
            return DelayedGetitem(obj=self, key=key)

    def __setitem__(self, key, value):
        if key in self.value:
            raise RuntimeError(f"store is immutable, cannot overwrite key {key}")
        if isinstance(value, DelayedCallable):
            self.queue.append(value)
        self.value[key] = value

    def __delitem__(self, key):
        raise RuntimeError(f"store is immutable, cannot delete key {key}")

    def __contains__(self, key):
        return key in self.value

    def __iter__(self):
        return iter(self.value)

    def clear(self):
        return self.value.clear()

    def __copy__(self):
        new = type(self).__new__()
        new.value = self.value.copy()
        return new

    copy = __copy__

    @classmethod
    def fromkeys(cls, iterable, value):
        raise NotImplementedError

    def get(self, key, default):
        if key in self:
            return self[key]
        else:
            return default

    def items(self):
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
    def __init__(self, queue, output_path):
        self.output_path = output_path
        super().__init__(queue)
