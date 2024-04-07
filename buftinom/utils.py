from typing import Any, Callable


def filter_max(iterable, key_func: Callable[[Any], float]):
    if len(iterable) == 0:
        return 0

    weighted = [(i, key_func(i)) for i in iterable]
    max_val = max([w for i, w in weighted])
    return [i for i, w in list(filter(lambda v: v[1] == max_val, weighted))]
