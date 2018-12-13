#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import os
import pickle
from io import BufferedReader, BufferedWriter
from gzip import GzipFile
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
            with BufferedReader(GzipFile(b["path"], "rb")) as fh:
                while True:
                    try:
                        k, v = pickle.load(fh)
                    except EOFError:
                        break
                    else:
                        yield k, v

    @property
    def size(self) -> int:
        size = 0
        for b in self.buckets:
            try:
                s = os.path.getsize(b["path"])
            except FileNotFoundError:
                continue
            else:
                size += s

        return size

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
                with BufferedWriter(GzipFile(b["path"], "ab")) as fh:
                    pickle.dump(b["data"], fh)
                b["data"] = {}

    @staticmethod
    def merge_bucket(bucket: dict, as_list: bool=False):
        data = {}
        if os.path.isfile(bucket["path"]):
            with BufferedReader(GzipFile(bucket["path"], "rb")) as fh:
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

    def merge(self, as_list: bool=False) -> Generator[dict, None, None]:
        for b in self.buckets:
            yield self.merge_bucket(b, as_list)

    def remove(self):
        for b in self.buckets:
            try:
                os.remove(b["path"])
            except FileNotFoundError:
                pass

        os.rmdir(self.path)


class Store(object):
    def __init__(self, path: Optional[str]=None, dir=None):
        if path:
            self.path = path
            if os.path.isfile(self.path):
                self.fh = BufferedReader(GzipFile(self.path, "rb"))
            else:
                self.fh = BufferedWriter(GzipFile(self.path, "wb"))
            self.delete_on_close = False
        else:
            fd, self.path = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.path)
            self.fh = BufferedWriter(GzipFile(self.path, "wb"))
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
        self.fh = BufferedReader(GzipFile(self.path, "rb"))
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
