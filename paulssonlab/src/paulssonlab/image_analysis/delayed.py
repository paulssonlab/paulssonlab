import itertools as it
from collections.abc import Callable, Hashable, MutableMapping, MutableSequence
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any

import dask
import distributed

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
    args: list
    kwargs: dict
    dependencies: "list[Delayed]" = field(default_factory=list)
    runner: "Runner"
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
            args = unbox_delayed(self.args)
            kwargs = unbox_delayed(self.kwargs)
            res = self.runner.run(self.func, *args, **kwargs)
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


class Runner:
    def __init__(self):
        self.queue = []

    @staticmethod
    def get(config):
        if config is False or config == "eager":
            return EagerRunner()
        elif config is True or config == "delayed":
            return DaskRunner()
        elif isinstance(config, distributed.Client):
            return DaskDistributedRunner(config)
        else:
            raise ValueError(
                f"unknown Runner config {config}, must be eager, delayed, or a distributed.Client instance"
            )

    def delayed(self, func, *args, **kwargs):
        dependencies = list(get_delayed(func, args, kwargs))
        if dependencies:
            delayed = DelayedCallable(
                func=func,
                args=args,
                kwargs=kwargs,
                dependencies=dependencies,
                runner=self,
            )
            self.queue.append(delayed)
            return delayed
        else:
            return self.run(func, *args, **kwargs)

    def run(self, func, *args, **kwargs):
        raise NotImplementedError

    def poll(self):
        while True:
            any_fired = False
            for idx in enumerate(self.queue):
                delayed = self.queue[idx]
                if delayed.is_ready():
                    delayed.result()
                    del self.queue[idx]
                    any_fired = True

            if not any_fired:
                break
        return


class EagerRunner(Runner):
    def run(self, func, *args, **kwargs):
        return func(*args, **kwargs)


class DaskRunner(Runner):
    def run(self, func, *args, **kwargs):
        return dask.delayed(func)(*args, **kwargs)


class DaskDistributedRunner(Runner):
    def __init__(self, client):
        super().__init__()
        self.client = client

    def run(self, func, *args, **kwargs):
        return self.client.submit(func, *args, **kwargs)


class DelayedStore:
    def __init__(self):
        self._store = {}

    def __len__(self):
        return len(self._store)

    def __getitem__(self, key):
        if key in self:
            return self._store[key]
        else:
            return DelayedGetitem(obj=self, key=key)

    def __setitem__(self, key, value):
        if key in self:
            raise RuntimeError(f"store is immutable, cannot overwrite key {key}")
        self._store[key] = value

    def __delitem__(self, key):
        raise RuntimeError(f"store is immutable, cannot delete key {key}")

    def __contains__(self, key):
        return key in self._store

    def __iter__(self):
        return iter(self._store)

    def clear(self):
        return self._store.clear()

    def __copy__(self):
        new = type(self).__new__()
        new._store = self._store.copy()
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
        return self._store.items()

    def keys(self):
        return self._store.keys()

    def pop(self, key, default=None):
        raise NotImplementedError

    def popitem(self):
        raise NotImplementedError

    def __reversed__(self):
        return reversed(self._queue)

    def setdefault(self, key, value):
        if key in self:
            return self[key]
        else:
            self[key] = value
            return value

    def update(self, *args, **kwargs):
        raise NotImplementedError

    def values(self):
        return self._store.values()

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
