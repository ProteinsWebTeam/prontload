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

    def merge(self, key):
        i = self.keys.index(key)
        b = self.buckets[i]

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

        return data
