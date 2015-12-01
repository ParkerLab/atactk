#
# atactk: ATAC-seq toolkit
#
# Copyright 2015 The Parker Lab at the University of Michigan
#
# Licensed under Version 3 of the GPL or any later version
#

"""
Utility code used in atactk.
"""

import collections
import operator


def add_lists(l1, l2):
    """
    Adds the values of two lists, entrywise.

    >>> add_lists([0, 1, 2], [3, 4, 5])
    [3, 5, 7]

    Parameters
    ----------
    l1: list
       The first list.
    l2: list
       The second list.

    Returns
    -------
    sum: list
        The list of the entrywise sums of the two lists' elements.

    """
    return map(operator.add, l1, l2)


def take(count, seq):
    """
    Return a list of up to `count` elements from the iterable `seq`.

    Parameters
    ----------
    count: int
        The number of elements to take from `seq`.
    seq: iterator-or-iterable
        An iterator or iterable from which to take elements.

    Returns
    -------
    list
        A list of up to `count` elements. There may be fewer if `seq` has been exhausted.

    """
    if not isinstance(seq, collections.Iterator):
        seq = iter(seq)
    l = []
    try:
        for i in range(count):
            l.append(next(seq))
    except StopIteration:
        pass
    return l


def partition(count, seq):
    """
    Create a generator of lists of `count` elements from `seq`.

    >>> list(partition(3, range(1, 10)))
    [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

    If `seq` isn't a multiple of `count`, the last list will contain
    the remaining items.

    >>> list(partition(3, range(1, 9)))
    [[1, 2, 3], [4, 5, 6], [7, 8]]

    Parameters
    ----------
    count: int
        The number of elements of `seq` to put in each partition.
    seq: iterator-or-iterable
        An iterator or iterable to be partitioned.

    Yields
    ------
    list
        A list representing a partition of `count` elements.
    """

    if not isinstance(seq, collections.Iterator):
        seq = iter(seq)
    partition = take(count, seq)
    while partition:
        yield partition
        partition = take(count, seq)
