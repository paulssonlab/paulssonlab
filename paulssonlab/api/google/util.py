def paginator(req_func):
    @wraps(req_func)
    def f(*args, **kwargs):
        res = req_func(*args, **kwargs).execute()
        for item in res["items"]:
            yield item
        while "nextPageToken" in res and res["nextPageToken"]:
            res = req_func(
                *args, **{"pageToken": res["nextPageToken"], **kwargs}
            ).execute()
            for item in res["items"]:
                yield item

    return f
