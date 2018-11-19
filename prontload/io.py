#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import os
import pickle
from tempfile import mkdtemp, mkstemp


class Organiser(object):
    def __init__(self, keys, path=None, dir=None):
        self.keys = keys

        if path is None:
            self.path = mkdtemp(dir=dir)
        else:
            self.path = path
            os.makedirs(self.path, exist_ok=True)

        self.buckets = [
            {
                "path": os.path.join(self.path, str(key)),
                "data": {}
            }
            for key in self.keys
        ]
        self.index = 0

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        try:
            b = self.buckets[self.index]
        except IndexError:
            raise StopIteration
        else:
            self.index += 1
            with open(b["path"], "rb") as fh:
                return pickle.load(fh)

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
                with open(b["path"], "ab") as fh:
                    pickle.dump(b["data"], fh)
                b["data"] = {}

    def merge(self):
        size_before = 0
        size_after = 0

        for b in self.buckets:
            size_before += os.path.getsize(b["path"])

            data = {}
            with open(b["path"], "rb") as fh:
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

            with open(b["path"], "wb") as fh:
                pickle.dump(data, fh)

            size_after += os.path.getsize(b["path"])

        return size_before, size_after

    def remove(self):
        for b in self.buckets:
            os.remove(b["path"])
        os.rmdir(self.path)


class Store(object):
    def __init__(self, path=None, dir=None, mode="rb"):
        if path:
            self.path = path
            self.temporary = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.temporary = True

        self.fh = open(self.path, mode)

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
        self.fh = open(self.path, "rb")
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
        # Expects `self.path` to be open for writing
        pickle.dump(obj, self.fh)

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None

    def remove(self):
        if self.temporary and self.path is not None:
            os.remove(self.path)
            self.path = None
