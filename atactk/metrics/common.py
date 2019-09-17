#
# atactk: ATAC-seq toolkit
#
# Copyright 2015 Stephen Parker
#
# Licensed under Version 3 of the GPL or any later version
#


import atactk.data
import atactk.util


DEFAULT_CUT_POINT_OFFSET = 4


def reduce_scores(scores, resolution):
    """
    Reduce a sequence of scores by summing every `resolution` values.

    Called with `scores` of [0, 1, 1, 4, 2], you'd get the following
    results at various resolutions:

    ==========  ======
    Resolution  Result
    ==========  ======
    1           [0, 1, 1, 4, 2]
    2           [1, 5, 2]
    3           [2, 6]
    4           [6, 2]
    10          [8]
    ==========  ======
    """
    if resolution == 1:
        return scores
    return [sum(chunk) for chunk in atactk.util.partition(resolution, scores)]
