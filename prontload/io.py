#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import gzip
import os
import pickle
from tempfile import mkdtemp, mkstemp
from typing import Union


class Organiser(object):
    def __init__(self, keys, path=None, dir=None, compress=True):
        self.keys = keys

        if path is None:
            self.path = mkdtemp(dir=dir)
        else:
            self.path = path
            os.makedirs(self.path, exist_ok=True)

        self.open = gzip.open if compress else open

        self.buckets = [
            {
                "path": os.path.join(self.path, str(i+1)),
                "data": {}
            }
            for i in range(len(self.keys))
        ]

    def __iter__(self):
        for b in self.buckets:
            with self.open(b["path"], "rb") as fh:
                while True:
                    try:
                        k, v = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        yield k, v

    @property
    def size(self):
        return len(self.buckets)

    def add(self, key, value):
        i = bisect.bisect_right(self.keys, key)
        if i:
            bucket = self.buckets[i-1]
            if key in bucket["data"]:
                bucket["data"][key].append(value)
            else:
                bucket["data"][key] = [value]
        else:
            raise ValueError(key)

    def dump(self):
        for b in self.buckets:
            if b["data"]:
                with self.open(b["path"], "ab") as fh:
                    pickle.dump(b["data"], fh)
                b["data"] = {}

    def merge(self):
        size_before = 0
        size_after = 0

        for b in self.buckets:
            data = {}
            if os.path.isfile(b["path"]):
                size_before += os.path.getsize(b["path"])
                with self.open(b["path"], "rb") as fh:
                    while True:
                        try:
                            chunk = pickle.load(fh)
                        except EOFError:
                            break
                        else:
                            for key, value in chunk.items():
                                if key in data:
                                    data[key] += value
                                else:
                                    data[key] = value

            with self.open(b["path"], "wb") as fh:
                for key in sorted(data):
                    pickle.dump((key, data[key]), fh)

            size_after += os.path.getsize(b["path"])

        return max(size_before, size_after)

    def remove(self):
        for b in self.buckets:
            os.remove(b["path"])
        os.rmdir(self.path)


class Store(object):
    def __init__(self, path: Union[str, None]=None, dir=None):
        if path:
            self.path = path
            self.fh = gzip.open(self.path, "rb")
            self.delete_on_close = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.fh = gzip.open(self.path, "wb")
            self.delete_on_close = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        self.remove()

    def __del__(self):
        self.close()
        self.remove()

    def __iter__(self):
        self.close()
        self.fh = gzip.open(self.path, "rb")
        return self

    def __next__(self):
        try:
            obj = pickle.load(self.fh)
        except EOFError:
            raise StopIteration
        else:
            return obj

    @property
    def size(self):
        try:
            return os.path.getsize(self.path)
        except (FileNotFoundError, TypeError):
            return 0

    def add(self, obj):
        pickle.dump(obj, self.fh)

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None

    def remove(self):
        if self.path is not None and self.delete_on_close:
            os.remove(self.path)
            self.path = None
