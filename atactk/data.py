#
# atactk: ATAC-seq toolkit
#
# Copyright 2015 Stephen Parker
#
# Licensed under Version 3 of the GPL or any later version
#


"""
Code for reading and manipulating data commonly used in ATAC-seq pipelines.
"""

import csv
import gzip
import pysam
import sys


FEATURE_FIELDNAMES = [
    'reference',
    'start',
    'end',
    'name',
    'score',
    'strand',
]

NUCLEOTIDE_COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "N": "N",
    "T": "A",
    "a": "t",
    "c": "g",
    "g": "c",
    "n": "n",
    "t": "a",
}


class ExtendedFeature(object):
    """
    A feature plus a fixed extended region.

    You can define the region by passing the `extension` parameter to the constructor, e.g.::

        feature = ExtendedFeature(extension=100, **bed_record)

    Most of :class:`ExtendedFeature`'s attributes map to the first six
    fields in a BED file. Where our names for the fields differ, the
    BED format name from https://genome.ucsc.edu/FAQ/FAQformat.html is
    included in parentheses below.

    Attributes
    ----------

    reference: str
        The reference sequence on which the feature is located. (``chrom``)
    feature_start: int
        The starting position of the feature in the reference sequence, zero-based. (``chromStart``)
    feature_end: int
        The ending position of the feature in the reference sequence, which is one past the last base in the feature. (``chromEnd``)
    name: str
        The name of the feature.
    score: float
        A numeric score.
    strand: str
        Either ``+`` or ``-``.
    """

    def __init__(self, reference=None, start=None, end=None, name=None, score=0, strand=None, extension=100):

        # required BED fields
        self.reference = reference
        self.feature_start = int(start)
        self.feature_end = int(end)

        # optional BED fields
        self.name = name or None
        self.score = score and float(score) or None
        self.strand = strand or None

        # region adjustments
        self.extension = int(extension)

    def __str__(self):
        return '\t'.join(str(attribute or '') for attribute in [
            self.reference,
            self.feature_start,
            self.feature_end,
            self.name,
            self.score,
            self.strand,
            self.extension,
        ])

    @property
    def center(self):
        return self.feature_start + int(round(self.feature_length / 2.0))

    @property
    def feature_length(self):
        return self.feature_end - self.feature_start

    @property
    def is_reverse(self):
        return self.strand == '-'

    @property
    def region_end(self):
        return self.center + self.extension

    @property
    def region_length(self):
        return self.extension * 2

    @property
    def region_start(self):
        return self.center - self.extension


def complement(seq):
    """
    Return the complement of the supplied nucleic sequence.

    Nucleic of course implies that the only recognized bases are A, C,
    G, T and N. Case will be preserved.

    Parameters
    ----------
    seq: str
        A nucleic sequence.

    Returns
    -------
    str
        The complement of the given sequence.
    """
    return ''.join(NUCLEOTIDE_COMPLEMENTS[base] for base in seq)


def reverse_complement(seq):
    """
    Return the reverse complement of the supplied nucleic sequence.

    Parameters
    ----------
    seq: str
        A nucleic sequence.

    Returns
    -------
    str
        The reverse complement of the given sequence.

    See also
    --------
    :func:`~atactk.data.complement`
    """
    return complement(reversed(seq))


def open_maybe_gzipped(filename):
    """
    Open a possibly gzipped file.

    Parameters
    ----------
    filename: str
        The name of the file to open.

    Returns
    -------
    file
        An open file object.
    """
    with open(filename, 'rb') as test_read:
        byte1, byte2 = ord(test_read.read(1)), ord(test_read.read(1))
        if byte1 == 0x1f and byte2 == 0x8b:
            f = gzip.open(filename, mode='rt')
        else:
            f = open(filename, 'rt')
    return f


def count_features(filename):
    count = 0
    with open_maybe_gzipped(filename) as f:
        for line in f:
            count += 1
    return count


def read_features(filename, extension=100, feature_class=ExtendedFeature):
    """
    Return a generator of :class:`ExtendedFeature` instances from the named tab-separated value file.

    Most BED-like files should work; we read the three required and
    first three optional BED fields to get coordinates, and any extra
    fields are ignored.

    Parameters
    ----------
    filename: str
        The (optionally gzipped) tab-separated value file from which to read features. Use '-' to read from standard input.
    extension: int
        The number of bases to score on either side of each feature.
    feature_class: class
        Each row of the file will be instantiated with this class.

    Yields
    ------
    feature
        An :class:`ExtendedFeature` instance for each row of the file.
    """

    with filename == '-' and sys.stdin or open_maybe_gzipped(filename) as source:
        reader = csv.DictReader(source, fieldnames=FEATURE_FIELDNAMES, restkey='extra_fields', dialect='excel-tab')

        for row in reader:
            if 'extra_fields' in row:
                del row['extra_fields']
            yield feature_class(extension=extension, **row)


def open_alignment_file(alignment_filename):
    alignment_file = pysam.AlignmentFile(alignment_filename, 'rb')
    try:
        alignment_file.check_index()
    except AttributeError:
        raise AttributeError('The alignments file {} is not in BAM format. Please supply an indexed BAM file.'.format(alignment_filename))
    except ValueError:
        raise ValueError('The alignment file {} is not usable. Please supply an indexed BAM file.'.format(alignment_filename))

    return alignment_file


def filter_aligned_segments(aligned_segments, include_flags, exclude_flags, quality):
    """
    Filter aligned segments using SAM flags and mapping quality.

    Parameters
    ----------
    aligned_segments: list
        Aligned reads to filter.
    include_flags: list
        Reads matching any include flag will be returned.
    exclude_flags: list
        Reads matching any exclude flag will not be returned.
    quality: int
        Only reads with at least this mapping quality will be returned.

    Returns
    -------
    filtered_aligned_segments: list
        The set of the aligned segments supplied to the function which
        meet the specified criteria.

    Examples
    --------

    You probably want `include_flags` of [83, 99, 147, 163] and
    `exclude_flags` of [4, 8].

    Flag 4 means the read is unmapped, 8 means the mate is unmapped.

    Properly paired and mapped forward aligned segments have flags in [99, 163]

    99:
       -   1: read paired
       -   2: read mapped in proper pair
       -  32: mate reverse strand
       -  64: first in pair

    163:
       -   1: read paired
       -   2: read mapped in proper pair
       -  32: mate reverse strand
       - 128: second in pair

    Properly paired and mapped reverse aligned segments have flags in [83, 147].

    83:
       -   1: read paired
       -   2: read mapped in proper pair
       -  16: read reverse strand
       -  64: first in pair

    147:
       -   1: read paired
       -   2: read mapped in proper pair
       -  16: read reverse strand
       - 128: second in pair
    """

    filtered_aligned_segments = [a for a in aligned_segments if all([
        a.mapping_quality >= quality,
        any(map(lambda f: (a.flag & f) == f, include_flags)),
        all(map(lambda f: (a.flag & f) == 0, exclude_flags))
    ])]
    return filtered_aligned_segments


def make_fastq_pair_reader(fastq_file1, fastq_file2):
    """
    Return a generator producing pairs of records from two FASTQ files.

    The intent is to produce read pairs from paired-end sequence data.

    Parameters
    ----------
    fastq_file1: str
        The name of the first FASTQ file.

    fastq_file2: str
        The name of the second FASTQ file.

    Yields
    ------
    tuple
        A tuple containing two 4-element lists, one for each FASTQ
        record, representing the ID, sequence, comment, and quality lines.
    """

    f1 = open_maybe_gzipped(fastq_file1)
    f2 = open_maybe_gzipped(fastq_file2)
    while True:
        yield (
            [
                next(f1).strip(),  # name
                next(f1).strip(),  # sequence
                next(f1).strip(),  # comment ('+' line)
                next(f1).strip()   # quality
            ],
            [
                next(f2).strip(),  # name
                next(f2).strip(),  # sequence
                next(f2).strip(),  # comment ('+' line)
                next(f2).strip()   # quality
            ],
        )
