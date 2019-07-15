import numpy as np
from numba import int64, jitclass


@jitclass(
    [
        ("max_size", int64),
        ("front", int64),
        ("back", int64),
        ("i", int64[:]),
        ("j", int64[:]),
    ]
)
class Fifo:
    def __init__(self, max_size):
        self.max_size = max_size
        self.front = 0
        self.back = 0
        self.i = np.empty((max_size,), dtype=np.int_)
        self.j = np.empty((max_size,), dtype=np.int_)

    def _wrap(self, idx):
        return idx % self.max_size

    @property
    def empty(self):
        return self.back == self.front

    def clear(self):
        self.front = 0
        self.back = 0

    def push(self, i, j):
        self.i[self._wrap(self.back)] = i
        self.j[self._wrap(self.back)] = j
        self.back += 1

    def pop(self):
        i = self.i[self._wrap(self.front)]
        j = self.j[self._wrap(self.front)]
        self.front += 1
        if self.empty:
            self.clear()
        return i, j
