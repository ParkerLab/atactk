#
# atactk: ATAC-seq toolkit
#
# Copyright 2015 The Parker Lab at the University of Michigan
#
# Licensed under Version 3 of the GPL or any later version
#


"""
Code for making quantitative observations about ATAC-seq experiments.
"""

import collections
import functools

import pysam

import atactk.data
import atactk.util


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


def aggregate_scores(scores, extension, resolution):
    """
    Adjust scores in the extended region around a feature.

    Given a sequence containing the score at each base in a region,
    the size of the extended region around the feature, and the
    desired resolution in that extended region, reduce the extended
    scores.

    Parameters
    ----------
    scores: list
        A list containing a score for each base in a region around a feature.
    extension: int
        The number of bases at the beginning and end of the list considered the extended region.
    resolution: int
        The desired scoring resolution in the extended region.

    See Also
    --------
    reduce_scores: Reduce scores by summing every `resolution` values.

    """
    return (
        reduce_scores(scores[:extension], resolution) +
        scores[extension:-extension] +
        reduce_scores(scores[-extension:], resolution)
    )


def find_cut_point(aligned_segment, cut_point_offset=4):
    """Return the position of the given aligned segment's ATAC-seq cut point.

    Parameters
    ----------
    aligned_segment: :class:`pysam.AlignedSegment`
        https://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment

    Returns
    -------
    int
        Position of the ATAC-seq cut point.
    """
    if aligned_segment.is_reverse:
        # the cut point is the reference_end minus (cut_point_offset + 1)
        # (pysam's reference_end is one past the last aligned residue)
        cut_point = aligned_segment.reference_end - (cut_point_offset + 1)
    else:
        cut_point = aligned_segment.reference_start + cut_point_offset  # start of the read plus offset
    return cut_point


def count_cut_points(aligned_segments, start, end, cut_point_offset=4):
    """
    Return any cut points in the region from the aligned segments.

    The cut point is the fifth base of the read, after the transposase
    integration. This implies of course that any aligned segments shorter than
    five bases cannot have a cut point and won't be part of the
    returned values.

    Parameters
    ----------
    aligned_segments: list
        A list of :class:`pysam.AlignedSegment`.
    start: int
        The start of the region of interest.
    end: int
        The end of the region of interest.

    Returns
    -------
    list
        A list of counts, one for each position from `start` to `end`, of cut points in the aligned segments that fell between the `start` and `end`..

    """

    cut_points_in_region = []
    for segment in aligned_segments:
        cut_point = find_cut_point(segment, cut_point_offset)
        if start <= cut_point < end:
            cut_points_in_region.append(cut_point)

    # initialize the region with zero counts
    counts = {p: 0 for p in range(start, end)}

    # add the cut points
    counts.update(collections.Counter(cut_points_in_region))

    cut_point_counts = [counts[position] for position in sorted(counts.keys())]
    return cut_point_counts


def add_cut_points_to_region_tree(region_tree, group_key, strand, cut_points):
    """
    Record cut points by position, template size, and strand.

    The tree consists of nested dictionaries. The first level keys are
    positions in a region. The second level keys are the names of
    template size groups. The third level keys are the strand, and the
    values are the cut point counts for that position, template size
    group, and strand.


    Parameters
    ----------
    region_tree: dict
        The tree to which the cut point counts will be added.
    group_key: str
        The key corresponding to the template size group.
    strand: str
        The strand of the cut point's aligned segment: ``+`` or ``-``.
    cut_points: list
        The count of cut points at each position in the region that matched the template size group and strand.

    Notes
    -----

    Only positive counts are recorded. It takes a lot of time and
    space to record so many zeroes, and it's better to produce them on
    demand via :class:`collections.defaultdict`. So instead, collect
    all the scores, and after that work is done, update a tree built
    with :class:`collections.defaultdict`, then work with that. See
    the ``make_cut_matrix`` script included with ``atactk`` for an
    example.
    """

    for position, count in enumerate(cut_points, 0 - int(len(cut_points) / 2)):
        if count > 0:
            if position not in region_tree:
                region_tree[position] = {}
            if group_key not in region_tree[position]:
                region_tree[position][group_key] = {}
            if strand not in region_tree[position][group_key]:
                region_tree[position][group_key][strand] = count
            else:
                region_tree[position][group_key][strand] += count


def score_feature(alignment_filename, bin_groups, include_flags, exclude_flags, quality, feature, cut_point_offset=4):
    """
    Count the number of transposition events around the given feature.

    Parameters
    ----------
    alignment_filename: str
        The BAM file containing aligned reads.
    bin_groups: iterable
        A sequence of iterables containing bins and the resolution with which they should be scored.
    include_flags: iterable
        The SAM flags to use when selecting aligned segments to score.
    exclude_flags: iterable
        The SAM flags to use when excluding aligned segments to score; any flag present on a read excludes it.
    quality: int
        The minimum mapping quality a read must have to be scored.
    feature: ExtendedFeature
        The feature to score.

    Returns
    -------
    tuple
        A tuple of `(row, tree)` where

        * `row` is a tab-separated list of scores in the region around the feature
        * `tree` is a three-level dict holding a score for each position, in each of the template size bins given, on each strand, e.g.::

            >>> tree[0]['36_149']['F']
            22
            >>> tree[0]['36_149']['R']
            15


    See Also
    --------
    add_cut_points_to_region_tree: Where the tree for the aggregate matrix is described more fully.

    """

    alignment_file = pysam.AlignmentFile(alignment_filename, 'rb')
    aligned_segments = alignment_file.fetch(feature.reference, max(0, feature.region_start), feature.region_end)
    aligned_segments = atactk.data.filter_aligned_segments(aligned_segments, include_flags, exclude_flags, quality)

    row = []
    tree = {}

    for group in bin_groups:
        group_rows = []
        group_key = '_'.join('%s.%s' % (bin[0], bin[1]) for bin in group)
        for (minimum_length, maximum_length, resolution) in group:
            bin_scores = []
            aligned_segments_in_bin = (a for a in aligned_segments if minimum_length <= abs(a.isize) <= maximum_length)
            forward_aligned_segments = (a for a in aligned_segments_in_bin if not a.is_reverse)
            reverse_aligned_segments = (a for a in aligned_segments_in_bin if a.is_reverse)

            forward_cut_points = count_cut_points(forward_aligned_segments, feature.region_start, feature.region_end, cut_point_offset)
            reverse_cut_points = count_cut_points(reverse_aligned_segments, feature.region_start, feature.region_end, cut_point_offset)

            if feature.is_reverse:
                forward_cut_points, reverse_cut_points = list(reversed(reverse_cut_points)), list(reversed(forward_cut_points))

            # for the discrete matrix: scores for each feature
            bin_scores.extend(aggregate_scores(forward_cut_points, feature.extension, resolution))
            bin_scores.extend(aggregate_scores(reverse_cut_points, feature.extension, resolution))
            group_rows.append(bin_scores)

            # for the aggregate matrix: scores for the entire region
            add_cut_points_to_region_tree(tree, group_key, 'F', forward_cut_points)
            add_cut_points_to_region_tree(tree, group_key, 'R', reverse_cut_points)

        if len(group) == 1:
            row.extend(group_rows[0])
        else:
            row.extend(functools.reduce(atactk.util.add_lists, group_rows))

    row = '\t'.join(str(score) for score in row)
    return row, tree
