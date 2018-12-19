#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import gzip
import os
import pickle
from tempfile import mkdtemp, mkstemp
from typing import Generator, List, Optional


class Organiser(object):
    def __init__(self, keys: List, dir: Optional[str]=None):
        self.keys = keys
        self.path = mkdtemp(dir=dir)
        self.buckets = [
            {
                "path": os.path.join(self.path, str(i+1)),
                "data": {}
            }
            for i in range(len(self.keys))
        ]

    def __iter__(self):
        for b in self.buckets:
            with gzip.open(b["path"], "rb") as fh:
                while True:
                    try:
                        k, v = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        yield k, v

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
                with gzip.open(b["path"], "ab") as fh:
                    pickle.dump(b["data"], fh)
                b["data"] = {}

    @staticmethod
    def merge_bucket(bucket: dict, as_list: bool=False):
        data = {}
        if os.path.isfile(bucket["path"]):
            with gzip.open(bucket["path"], "rb") as fh:
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

        if as_list:
            return [(key, data[key]) for key in sorted(data)]
        else:
            return data

    def merge(self) -> int:
        size_before = 0
        sife_after = 0
        for b in self.buckets:
            data = {}
            if os.path.isfile(b["path"]):
                size_before += os.path.getsize(b["path"])
                with gzip.open(b["path"], "rb") as fh:
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

            with gzip.open(b["path"], "wb") as fh:
                for key in sorted(data):
                    pickle.dump((key, data[key]), fh)

            sife_after += os.path.getsize(b["path"])

        return max(size_before, sife_after)

    def remove(self):
        for b in self.buckets:
            try:
                os.remove(b["path"])
            except FileNotFoundError:
                pass

        os.rmdir(self.path)


class Store(object):
    def __init__(self, path: Optional[str]=None, write: bool=False,
                 compresslevel: int=9, dir: Optional[str]=None):
        if path:
            self.path = path
            self.temporary = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            self.temporary = True
            write = True

        if write:
            self.fh = gzip.open(self.path, "wb", compresslevel)
        else:
            self.fh = gzip.open(self.path, "rb")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        if self.temporary:
            self.remove()

    def __del__(self):
        self.close()
        if self.temporary:
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
        except FileNotFoundError:
            return 0

    def add(self, obj):
        pickle.dump(obj, self.fh)

    def writefile(self, path: str):
        with open(path, "rb") as fh:
            self.fh.write(fh.read())
        os.remove(path)

    def close(self):
        self.fh.close()

    def remove(self):
        try:
            os.remove(self.path)
        except FileNotFoundError:
            pass
