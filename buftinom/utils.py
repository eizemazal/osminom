from typing import Any, Callable


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
    return [i for i, w in list(filter(lambda v: v[1] == max_val, weighted))]
