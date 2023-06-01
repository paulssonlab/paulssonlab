import json
import pickle
import base64
import zarr
import array
from collections.abc import Mapping
from functools import partial
import playhouse.sqlite_ext

# FROM: https://gist.github.com/simonw/7000493
class JSONEncoder(json.JSONEncoder):
    pickle_types = (zarr.Array, array.array)

    def default(self, obj):
        print(">>>", obj)
        # if isinstance(obj, Mapping):
        #     return { self.default(k): self.default(v) for k, v in obj.items() }
        if isinstance(obj, bytes):
            return {"__type__": "bytes", "value": base64.b64encode(obj).decode()}
        elif isinstance(obj, self.pickle_types):
            return {
                "__type__": "pickle",
                "value": base64.b64encode(pickle.dumps(obj)).decode(),
            }
        return super(JSONEncoder, self).default(obj)


class JSONDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, obj):
        print(">>>", obj)
        if "__type__" not in obj:
            return obj
        type = obj["__type__"]
        s = obj["value"]
        if type == "bytes":
            return base64.b64decode(s.encode())
        elif type == "pickle":
            return pickle.loads(base64.b64decode(s.encode()))
        return obj


dumps = partial(json.dumps, cls=JSONEncoder)
loads = partial(json.loads, cls=JSONDecoder)
JSONField = partial(playhouse.sqlite_ext.JSONField, json_dumps=dumps, json_loads=loads)

# class JSONField(playhouse.sqlite_ext.JSONField):
#     def python_value(self, value):
#         if value is not None:
#             try:
#                 return json.loads(value, cls=JSONDecoder)
#             except (TypeError, ValueError):
#                 return value

#     def db_value(self, value):
#         if value is not None:
#             return json.dumps(value, cls=JSONEncoder)


def demo():
    data = {"foo": {"bar": zarr.zeros(10)}}
    s = json.dumps(data, cls=JSONEncoder, indent=2)
    print(data)
    print(s)
    print(json.loads(s, cls=JSONDecoder))


if __name__ == "__main__":
    demo()
