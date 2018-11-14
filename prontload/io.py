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

    def load(self, key):
        i = self.keys.index(key)
        b = self.buckets[i]
        with open(b["path"], "rb") as fh:
            return pickle.load(fh)


class ProteinStore(object):
    def __init__(self, path=None, dir=None, mode="rb"):
        if path:
            self.path = path
            self.delete = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.delete = True

        self.fh = open(self.path, mode)

    def add(self, obj):
        # Expects `self.path` to be open for writing
        pickle.dump(obj, self.fh)

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None

        if self.delete and self.path is not None:
            os.remove(self.path)
            self.path = None

    def __del__(self):
        self.close()

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
