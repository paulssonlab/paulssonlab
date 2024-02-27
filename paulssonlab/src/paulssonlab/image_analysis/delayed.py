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


NoResult = object()


@dataclass(kw_only=True)
class DelayedValue(Delayed):
    value: Any = NoResult

    def is_ready(self):
        return self.value is not NoResult

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
    _result: Any = NoResult

    @cached_property
    def dependencies(self):
        # TODO: make sure this doesn't get called if dependencies is set
        return get_delayed(self.func, self.args, self.kwargs)

    def is_ready(self):
        return all(dep.is_ready() for dep in self.dependencies)

    @property
    def is_pending(self):
        return self._result is NoResult

    def result(self):
        if (res := self._result) is NoResult:
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
        return dask.delayed(func, *args, **kwargs)


class DaskDistributedRunner(Runner):
    def __init__(self, client):
        self.client = client

    def run(self, func, *args, **kwargs):
        return self.client.submit(func, *args, **kwargs)


class DelayedStore:
    def __init__(self, output_path):
        self.output_path = output_path
        self._store = {}

    def __contains__(self, key):
        return key in self._store

    def __getitem__(self, key):
        if key in self:
            return self._store[key]
        else:
            return DelayedGetitem(self, key)

    def __setitem__(self, key, value):
        if key in self:
            raise RuntimeError(f"attempting to overwrite key {value}")

    def getdefault(self, key, default):
        pass

    def setdefault(self, key, value):
        pass

    def write(self):
        raise NotImplementedError


class DelayedArrayStore(DelayedStore):
    pass


class DelayedTableStore(DelayedStore):
    pass
