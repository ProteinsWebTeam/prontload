#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import os
import pickle
from tempfile import mkdtemp


class Organiser(object):
    def __init__(self, keys, path=None, tmpdir=None):
        self.keys = keys

        if path is None:
            self.path = mkdtemp(dir=tmpdir)
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


class ProteinIterator(object):
    def __init__(self, filepath, mode="rb"):
        self.filepath = filepath
        self.mode = mode
        self.fh = None

    def add(self, obj):
        pickle.dump(obj, self.fh)

    def open(self, mode="rb"):
        self.close()
        self.fh = open(self.filepath, mode)

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None

    def __enter__(self):
        self.open(self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        try:
            obj = pickle.load(self.fh)
        except EOFError:
            raise StopIteration
        else:
            return obj
