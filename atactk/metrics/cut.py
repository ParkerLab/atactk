#
# atactk: ATAC-seq toolkit
#
# Copyright 2015 Stephen Parker
#
# Licensed under Version 3 of the GPL or any later version
#


"""
Code for making quantitative observations about ATAC-seq experiments.
"""

import collections
import functools
import multiprocessing
import signal

import atactk.data
import atactk.metrics.common as common
import atactk.util


def find_cut_point(aligned_segment, cut_point_offset=common.DEFAULT_CUT_POINT_OFFSET):
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


def count_cut_points(aligned_segments, start, end, cut_point_offset=common.DEFAULT_CUT_POINT_OFFSET):
    """
    Return any cut points in the region from the aligned segments.

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
        A list of counts, one for each position from `start` to `end`,
        of the aligned segments' cut points at that position.
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
    Record cut points by position, fragment size, and strand.

    The tree consists of nested dictionaries. The first level keys are
    positions in a region. The second level keys are the names of
    fragment size groups. The third level keys are the strand, and the
    values are the cut point counts for that position, fragment size
    group, and strand.


    Parameters
    ----------
    region_tree: dict
        The tree to which the cut point counts will be added.
    group_key: str
        The key corresponding to the fragment size group.
    strand: str
        The strand of the cut points' aligned segments: ``+`` or ``-``.
    cut_points: list
        The count of cut points at each position in the region that matched the fragment size group and strand.

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

    for position, count in enumerate(cut_points, 0 - (len(cut_points) // 2)):
        if count > 0:
            if position not in region_tree:
                region_tree[position] = {}
            if group_key not in region_tree[position]:
                region_tree[position][group_key] = {}
            if strand not in region_tree[position][group_key]:
                region_tree[position][group_key][strand] = count
            else:
                region_tree[position][group_key][strand] += count


def score_feature_with_cut_points(bin_groups, include_flags, exclude_flags, quality, cut_point_offset, feature):
    """
    Count the number of transposition events around the given feature.

    Intended to be run only within a `multiprocessing.Pool`, in which
    each worker initializes its own copy of the indexed BAM file
    containing aligned reads, stored in the global
    `atactk.metrics.cut.alignment_file`.

    Parameters
    ----------
    bin_groups: iterable
        A sequence of iterables containing fragment size bins and the resolution with which they should be scored. If
        omitted, the matrix will contain a column for each position in the extended region, representing the count of
        cuts at that position.
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
        * `tree` is a three-level dict holding a score for each position, in each of the fragment size bins given, on each strand, e.g.::

            >>> tree[0]['36_149']['F']
            22
            >>> tree[0]['36_149']['R']
            15


    See Also
    --------
    add_cut_points_to_region_tree: Where the tree for the cut point aggregate matrix is described more fully.
    """

    aligned_segments = alignment_file.fetch(feature.reference, max(0, feature.region_start), feature.region_end)
    aligned_segments = atactk.data.filter_aligned_segments(aligned_segments, include_flags, exclude_flags, quality)

    row = []
    tree = {}

    if bin_groups:
        for group in bin_groups:
            group_rows = []
            group_key = atactk.util.make_bin_group_key(group)
            for (minimum_length, maximum_length, resolution) in group:
                bin_scores = []
                aligned_segments_in_bin = [a for a in aligned_segments if minimum_length <= abs(a.isize) <= maximum_length]
                forward_aligned_segments = [a for a in aligned_segments_in_bin if not a.is_reverse]
                reverse_aligned_segments = [a for a in aligned_segments_in_bin if a.is_reverse]

                forward_cut_points = count_cut_points(forward_aligned_segments, feature.region_start, feature.region_end, cut_point_offset)
                reverse_cut_points = count_cut_points(reverse_aligned_segments, feature.region_start, feature.region_end, cut_point_offset)

                if feature.is_reverse:
                    # need to orient the cut point positions to the motif in the matrix
                    forward_cut_points, reverse_cut_points = list(reversed(reverse_cut_points)), list(reversed(forward_cut_points))

                # for the discrete matrix: scores for each feature
                bin_scores.extend(common.reduce_scores(forward_cut_points, resolution))
                bin_scores.extend(common.reduce_scores(reverse_cut_points, resolution))
                group_rows.append(bin_scores)

                # for the aggregate matrix: scores for the entire region
                add_cut_points_to_region_tree(tree, group_key, 'F', forward_cut_points)
                add_cut_points_to_region_tree(tree, group_key, 'R', reverse_cut_points)

            if len(group) == 1:
                row.extend(group_rows[0])
            else:
                row.extend(functools.reduce(atactk.util.add_lists, group_rows))
    else:
        row = count_cut_points(aligned_segments, feature.region_start, feature.region_end, cut_point_offset)
        if feature.is_reverse:
            row = list(reversed(row))
        add_cut_points_to_region_tree(tree, 'All', 'Both', row)

    row = '\t'.join(str(score) for score in row)
    return feature, row, tree


alignment_file = None


def scoring_process_init(alignment_filename):
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    global alignment_file
    if alignment_file is not None:
        alignment_file.close()
    alignment_file = atactk.data.open_alignment_file(alignment_filename)


def score_features_with_cut_points(alignment_file, bin_groups, include_flags, exclude_flags, quality, cut_point_offset, features, parallel=1):

    partial_scoring_function = functools.partial(
        score_feature_with_cut_points,
        bin_groups,
        include_flags,
        exclude_flags,
        quality,
        cut_point_offset
    )

    pool = multiprocessing.Pool(processes=parallel, initializer=scoring_process_init, initargs=[alignment_file])

    return pool.imap(partial_scoring_function, features, parallel)
