.. :changelog:

History
=======

0.1.4 (2015-12-17)
------------------

When generating an aggregate matrix, ensure that there is always a
line for each position, fragment size bin, and strand, even if there
is no signal there.

Support reading motifs from standard input, which required removing
time estimates, and returning the motif from
atactk.metrics.score_feature.

Remove option for reverse feature shift.


0.1.3 (2015-12-10)
------------------

Speed up scoring.

Open the alignment file once in each worker process, instead of in each
call to score_feature. There's a surprising amount of overhead in
pysam's opening of BAM files. The actual fetch calls are really quick.

The AlignmentFile instances are not supposed to be safe to share between
processes, so each worker still has to have its own, but it's a big
reduction in overhead; things are roughly twice as fast now.

Also refactor make_cut_matrix, mainly to enable profiling.

0.1.2 (2015-12-06)
------------------

Fix an overlooked realization of the list of motifs, which would be
copied into all processes. Use a generator expression for the results
when using only one process.

Convert a couple of other list comprehensions to generator expressions.

Change --parallel default to 1.

Improve logging.

0.1.1 (2015-12-01)
------------------

Arbitrary cut point location.

Now you can specify where (or whether) the cut point happens relative to
the ends of reads. We're using this to compare our footprinting with
DNase-seq data.

0.1.0 (2015-10-31)
------------------

* First release.
