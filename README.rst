.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4321161.svg
   :target: https://doi.org/10.5281/zenodo.4321161
======
atactk
======
A toolkit for working with ATAC-seq data.

What's in the box?
==================

Programs we've found useful in analyzing ATAC-seq data
------------------------------------------------------

* ``make_cut_matrix``: useful in conjunction with CENTIPEDE, and in
  generating plots of transcription factor binding sites.
* ``make_midpoint_matrix``: useful in conjunction with CENTIPEDE, and in
  generating plots of transcription factor binding sites.
* ``measure_signal``: measures size and position of ATAC-seq fragments
  overlapping genomic features.
* ``measure_features``: given a bigWig file of coverage counts and a BED
  file of features, calculates a requested statistic for each feature.
* ``plot_aggregate_cut_matrix.R``: generates plots for motifs given the
  aggregate output of `make_cut_matrix`
* ``plot_aggregate_midpoint_matrix.R``: generates plots for motifs given the
  aggregate output of `make_midpoint_matrix`
* ``plot_signal.R``: plots the output of `measure_signal`
* ``trim_adapters``: based on Jason Buenrostro's utility for trimming
  Illumina adapters by aligning paired reads to each other.

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

* Python. We've run it successfully under versions 2.7 and 3.5.
* pyBigWig
* pysam
* python-levenshtein
* sexpdata

Installation
============

At the command line::

  git clone https://github.com/ParkerLab/atactk
  pip install ./atactk

Or in one step, if you don't want a local copy::

  pip install git+https://github.com/ParkerLab/atactk

Documentation
=============

http://atactk.readthedocs.org/en/latest/

License
=======

GPLv3 or any later version.
