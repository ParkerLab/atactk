======
atactk
======

A toolkit for working with ATAC-seq data.

What's in the box?
==================

Programs we've found useful in ATAC-seq pipelines
-------------------------------------------------

* ``trim_adapters``: based on Jason Buenrostro's utility for trimming
  Illumina adapters by aligning paired reads to each other.
* ``make_cut_matrix``: useful in conjunction with CENTIPEDE, and in
  generating plots of transcription factor binding sites.
* ``plot_aggregate_matrix.R``: generates plots for motifs given the
  aggregate output of `make_cut_matrix`

A Python library you can use in your own tools for processing ATAC-seq data
---------------------------------------------------------------------------

The code underpinning our command-line tools has allowed us to make
our pipelines shorter and faster. Our ATAC-seq scoring functions work
directly with a BAM file and run in parallel, without the overhead of
invoking external applications. Particularly if you're trying to
produce quantitative metrics from your data, starting with your BAM
files, converting them to BED and bigWig so you can run bigWigSummary,
you might find your pipeline can be simplified too.

Requirements
============

* Python. We've run it successfully under versions 2.7.10 and 3.4.3.
* pysam
* python-levenshtein
* sexpdata

Installation
============

At the command line::

  git clone https://github.com/ParkerLab/atactk
  pip install ./atactk

Documentation
=============

http://atactk.readthedocs.org/en/latest/

License
=======

GPLv3 or any later version.
