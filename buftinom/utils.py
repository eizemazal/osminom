from typing import Any, Callable


def identity(x):
    return x


def filter_max(iterable, key_func: Callable[[Any], float]):
    """
    Returns items that score the maximum value in the key function.

    >>> filter_max([1, 4, 2, 4, 1], lambda x: x)
    [4, 4]
    >>> filter_max([-4, 1, 0, 2, 4, -3], lambda x: abs(x))
    [-4, 4]
    """
    if len(iterable) == 0:
        return []
    weighted = [(i, key_func(i)) for i in iterable]
    max_val = max([w for i, w in weighted])
    return [i for i, _ in list(filter(lambda v: v[1] == max_val, weighted))]


def choose(*, a, b, cmp: int):
    """Utility for choosing a or b or skip selection if cmp==0

    suggested use:

    selector = partial(choose, a='A', b='B')

    if result := selector(condition):
        return result

    if result := selector(condition2):
        return result
    """
    if cmp > 0:
        return a
    if cmp < 0:
        return b

    return None


def firstnn_idx(iterable, key: Callable[[Any], int]):
    for i, e in enumerate(iterable):
        if key(e) != 0:
            return i


def nonzero_indexes(iterable, key: Callable[[Any], int]):
    for i, e in enumerate(iterable):
        if key(e) != 0:
            yield i


def select(indexable, indexes):
    res = []
    for i in indexes:
        res.append(indexable[i])


def argmax_filter(iterable, include):
    m = 0
    for i, e in enumerate(iterable):
        if i not in include:
            continue

        if e != 0 and e > m:
            m = e

    return [i for i, e in enumerate(iterable) if e == m], m


def first_max(iterables: list[list[int]]):
    """
    See Iupac._fpod_many
    """
    row_indexes = set(range(len(iterables)))

    for column in zip(*iterables):
        rows, max_value = argmax_filter(column, row_indexes)
        if max_value == 0:
            continue

        rows = [row for row in rows if row in row_indexes]

        if len(rows) == 1:
            maxrow = rows[0]
            return maxrow

        row_indexes = set(rows)
