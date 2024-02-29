import itertools as it
from collections.abc import Callable, Hashable, MutableMapping, MutableSequence
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any

import dask
from dask.distributed import Client

from paulssonlab.util.core import iter_recursive, map_recursive


class DelayedNotReadyError(RuntimeError):
    pass


class Delayed:
    pass


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
    func: Callable
    args: list | None = None
    kwargs: dict | None = None
    dependencies: "list[Delayed]" = field(default_factory=list)
    _result: Any = NotReady

    @cached_property
    def dependencies(self):
        # TODO: make sure this doesn't get called if dependencies is set
        return get_delayed(self.func, self.args, self.kwargs)

    def is_ready(self):
        return all(dep.is_ready() for dep in self.dependencies)

    @property
    def is_pending(self):
        return self._result is NotReady

    def result(self):
        if (res := self._result) is NotReady:
            args = unbox_delayed(self.args) or ()
            kwargs = unbox_delayed(self.kwargs) or {}
            res = self.func(*args, **kwargs)
            self._result = res
        return res


@dataclass(kw_only=True)
class DelayedGetitem(Delayed):
    obj: MutableMapping | MutableSequence
    key: Hashable | int

    @property
    def dependencies(self):
        return get_delayed(self.obj, self.key, recursive=False)

    def is_ready(self):
        return self.key in self.obj

    @property
    def value(self):
        return self.obj[self.key]

    @value.setter
    def value(self, new_value):
        self.obj[self.key] = new_value

    def result(self):
        if not self.is_ready():
            raise DelayedNotReadyError
        return self.obj[self.key]


def get_delayed(*args, recursive=True):
    # TODO: look how dask implements unboxing
    if recursive:
        args = it.chain.from_iterable((iter_recursive(arg) for arg in args))
    return filter(lambda x: isinstance(x, Delayed), args)


def unbox_delayed(obj):
    return map_recursive(
        lambda x: x.result() if isinstance(x, Delayed) else x,
        obj,
    )


class DelayedQueue:
    def __init__(self):
        self._queue = []

    def delayed(self, func, *args, **kwargs):
        dependencies = list(get_delayed(func, args, kwargs))
        if dependencies:
            dc = DelayedCallable(
                func=func,
                args=args,
                kwargs=kwargs,
                dependencies=dependencies,
            )
            self.append(dc)
            return dc
        else:
            return func(*args, **kwargs)

    def poll(self):
        while True:
            any_fired = False
            idx = 0
            while idx < len(self._queue):
                delayed = self._queue[idx]
                if delayed.is_ready():
                    delayed.result()
                    del self._queue[idx]
                    any_fired = True
                else:
                    idx += 1

            if not any_fired:
                break
        return

    def append(self, delayed):
        self._queue.append(delayed)


class DelayedStore:
    def __init__(self, queue):
        self.value = {}
        self.queue = queue

    def __len__(self):
        return len(self.value)

    def __getitem__(self, key):
        if key in self:
            return self.value[key]
        else:
            return DelayedGetitem(obj=self, key=key)

    def __setitem__(self, key, value):
        if key in self:
            raise RuntimeError(f"store is immutable, cannot overwrite key {key}")
        if isinstance(value, Delayed):
            setter = DelayedCallable(
                func=self.value.__setitem__,
                args=[key, value],
                dependencies=[value],
            )
            self.queue.append(value)
            self.queue.append(setter)
        else:
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
        if key in self:
            return self[key]
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


class DelayedArrayStore(DelayedStore):
    pass


class DelayedTableStore(DelayedStore):
    pass
