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


#
# When determining a fragment's position relative to a feature, we consider it 'left' if both cut point positions are less than
# the feature's center reference position, 'right' if both are greater than the feature's center, and and 'overlapping' if they
# span the center position.
#
FEATURE_RELATIVE_POSITIONS = ('L', 'O', 'R')


def find_midpoint(aligned_segment):
    """
    Return the location of the midpoint of the aligned segment.

    For forward reads, it's reference_start plus half the fragment length.

    For reverse reads, it's reference_end plus half the fragment length.

    (Pysam stores fragment length as a negative number for reverse reads.)

    Parameters
    ----------
    aligned_segment: :class:`pysam.AlignedSegment`.

    Returns
    -------
    midpoint
        The integer position of the segment's midpoint.
    """
    if aligned_segment.is_reverse:
        return aligned_segment.reference_end + (aligned_segment.isize // 2)
    else:
        return aligned_segment.reference_start + (aligned_segment.isize // 2)


def count_midpoints(aligned_segments, feature):
    """
    Count the aligned segments' midpoints that fall in the feature's extended region.

    The list of counts at positions will be ordered with respect to
    the feature's strand. This just means that the returned sequence
    will be reversed for features on the reverse strand.

    Parameters
    ----------
    aligned_segments: list
        A list of :class:`pysam.AlignedSegment`.
    feature: ExtendedFeature
        The feature around which to count fragment midpoints.

    Returns
    -------
    generator
        A generator of counts, one for each position from `start` to `end`,
        of fragment midpoints at that position relative to the feature.
    """

    reads_seen = {}  # map of read names already seen, to avoid double-counting fragments
    midpoints_in_region = []
    for segment in aligned_segments:
        if segment.qname in reads_seen:
            # skip this one; we've already found the midpoint of its fragment using its mate
            continue

        reads_seen[segment.qname] = 1

        midpoint = find_midpoint(segment)
        if feature.region_start <= midpoint < feature.region_end:
            midpoints_in_region.append(midpoint)

    # initialize the region with zero counts
    counts = {p: 0 for p in range(feature.region_start, feature.region_end)}

    # add the cut points
    counts.update(collections.Counter(midpoints_in_region))

    midpoint_counts = [counts[position] for position in sorted(counts.keys())]
    if feature.is_reverse:
        midpoint_counts = list(reversed(midpoint_counts))
    return midpoint_counts


def add_midpoints_to_region_tree(region_tree, group_key, midpoints):
    """
    Record fragment midpoints by position and fragment size.

    The tree consists of nested dictionaries. The first level keys are
    positions in a region. The second level keys are the names of
    fragment size groups, and the values are the cut point counts for
    that position and fragment size group.


    Parameters
    ----------
    region_tree: dict
        The tree to which the cut point counts will be added.
    group_key: str
        The key corresponding to the fragment size group.
    midpoints: list
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

    for position, count in enumerate(midpoints, 0 - (len(midpoints) // 2)):
        if count > 0:
            if position not in region_tree:
                region_tree[position] = {}
            if group_key not in region_tree[position]:
                region_tree[position][group_key] = count
            else:
                region_tree[position][group_key] += count


def score_feature_with_midpoints(bin_groups, include_flags, exclude_flags, quality, feature):
    """
    Count the number of fragment midpoints around the given feature.

    Intended to be run only within a `multiprocessing.Pool`, in which
    each worker initializes its own copy of the indexed BAM file
    containing aligned reads, stored in the global
    `atactk.metrics.cut.alignment_file`.

    Parameters
    ----------
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
        * `tree` is a two-level dict holding a score for each position, in each of the fragment size bins given, e.g.::

            >>> tree[0]['36_149']
            22
            >>> tree[0]['36_149']
            15


    See Also
    --------
    add_midpoints_to_region_tree: Where the tree for the midpoint aggregate matrix is described more fully.

    """

    alignment_search_region_extension = max([b[1] for bins in bin_groups for b in bins]) // 2

    aligned_segments = alignment_file.fetch(
        feature.reference,
        max(0, feature.region_start - alignment_search_region_extension),
        feature.region_end + alignment_search_region_extension
    )
    aligned_segments = atactk.data.filter_aligned_segments(aligned_segments, include_flags, exclude_flags, quality)

    row = []
    tree = {}

    for group in bin_groups:
        group_rows = []
        group_key = atactk.util.make_bin_group_key(group)
        for (minimum_length, maximum_length, resolution) in group:
            bin_scores = []
            aligned_segments_in_bin = [a for a in aligned_segments if minimum_length <= abs(a.isize) <= maximum_length]
            midpoint_counts = count_midpoints(aligned_segments_in_bin, feature)

            # for the discrete matrix: scores for each feature
            bin_scores.extend(common.reduce_scores(midpoint_counts, resolution))
            group_rows.append(bin_scores)

            # for the aggregate matrix: scores for the entire region
            add_midpoints_to_region_tree(tree, group_key, midpoint_counts)

        if len(group) == 1:
            row.extend(group_rows[0])
        else:
            row.extend(functools.reduce(atactk.util.add_lists, group_rows))

    row = '\t'.join(str(score) for score in row)
    return feature, row, tree


alignment_file = None


def scoring_process_init(alignment_filename):
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    global alignment_file
    if alignment_file is not None:
        alignment_file.close()
    alignment_file = atactk.data.open_alignment_file(alignment_filename)


def score_features_with_midpoints(alignment_file, bin_groups, include_flags, exclude_flags, quality, cut_point_offset, features, parallel=1):

    partial_scoring_function = functools.partial(
        score_feature_with_midpoints,
        bin_groups,
        include_flags,
        exclude_flags,
        quality,
        cut_point_offset
    )

    pool = multiprocessing.Pool(processes=parallel, initializer=scoring_process_init, initargs=[alignment_file])

    return pool.imap(partial_scoring_function, features, parallel)


def get_fragment_cut_points(aligned_segment, cut_point_offset=common.DEFAULT_CUT_POINT_OFFSET):
    """
    Return the positions of both of the cut points in the the given aligned segment's fragment.

    80bp fragment mapped to reference:
    S-------------------------------------------------------------------------------E

    50bp paired end read 1 (forward):
    S---C---------------------------------------------E

    50bp paired end read 2 (reverse):
                                  S---------------------------------------------C---E

    Starts (S) and ends (E) are always just the AlignedSegment's reference_start and reference_end. The cut points are always
    relative to the start of the read, which for reverse reads is reference_end. Pysam's reference_end is one *past* the last
    actual position of the read, so we have to take that into account when calculating the cut point position.

    Parameters
    ----------
    aligned_segment: :class:`pysam.AlignedSegment`
        https://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment

    Returns
    -------
    tuple
        A tuple of (int, int) containing the fragment's sorted cut point positions.
    """
    if aligned_segment.is_reverse:
        right_cut_point = aligned_segment.reference_end - (cut_point_offset + 1)
        left_cut_point = aligned_segment.next_reference_start + cut_point_offset
    else:
        left_cut_point = aligned_segment.reference_start + cut_point_offset
        right_cut_point = aligned_segment.reference_start + aligned_segment.isize - (cut_point_offset + 1)
    return (left_cut_point, right_cut_point)


def get_fragment_position_relative_to_feature(aligned_segment, cut_point_offset, feature):
    """
    Determine where an aligned segment's fragment maps relative to a feature.

    If both cut points in the fragment map before the center position of the feature, we say the fragment is to the feature's
    left. If both map after the center, it's to the right. Otherwise the fragment overlaps the center.
    """

    cut_points = get_fragment_cut_points(aligned_segment, cut_point_offset)
    if feature.center > max(cut_points):
        # 'left' of the feature
        if feature.is_reverse:
            position = 'R'
        else:
            position = 'L'
    elif feature.center < min(cut_points):
        # 'right' of the feature
        if feature.is_reverse:
            position = 'L'
        else:
            position = 'R'
    else:
        # overlapping the feature
        position = 'O'
    return position


def find_midpoints_around_feature(alignment_filename, include_flags, exclude_flags, quality, cut_point_offset, feature):
    """
    Find fragment midpoints around a feature.

    Parameters
    ----------
    alignment_filename: str
        The BAM file containing aligned reads.
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
    A list containing for each midpoint found its position relative to the feature center, its fragment size, and its
    fragment's position relative to the feature center.
    """

    alignment_file = atactk.data.open_alignment_file(alignment_filename)

    feature_center = feature.center
    reads_seen = {}  # map of read names already seen, to avoid double-counting fragments
    midpoints = []

    aligned_segments = alignment_file.fetch(feature.reference, max(0, feature.region_start), feature.region_end)
    aligned_segments = atactk.data.filter_aligned_segments(aligned_segments, include_flags, exclude_flags, quality)

    for aligned_segment in aligned_segments:
        if aligned_segment.qname in reads_seen:
            # skip this one; we've already found the midpoint of its fragment using its mate
            continue

        reads_seen[aligned_segment.qname] = 1

        midpoint = find_midpoint(aligned_segment)
        if feature.region_start <= midpoint <= feature.region_end:
            distance_to_center = midpoint - feature_center
            if feature.is_reverse:
                distance_to_center *= -1
            fragment_relative_position = get_fragment_position_relative_to_feature(aligned_segment, cut_point_offset, feature)
            midpoints.append((feature_center, midpoint, distance_to_center, abs(aligned_segment.isize), fragment_relative_position))

    return midpoints
