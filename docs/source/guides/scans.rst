Included Scans
==============

The complete API reference for the included scans can be found at :doc:`bsmart.scans <../api/bsmart.scans>`.

They include:

* :doc:`Grid <../api/generated/bsmart.scans.Grid>`, a basic grid scan over a choice of parameter ranges. A very simple choice ideal for initial exploration of a couple of parameters, and very useful for collider studies. With this in mind, it can now be run using MPI, which allows points to be distributed across several nodes, each of which can then e.g. run MadGraph using several cores.
* :doc:`Random <../api/generated/bsmart.scans.Random>`, a basic random scan. A very simple choice for initial exploration of a few parameters.
* :doc:`MCMC <../api/generated/bsmart.scans.MCMC>`, a simple (and multicore) Markov-Chain Monte Carlo scan. This is the go-to scan for initial tests of a model.
* :doc:`read_csv <../api/generated/bsmart.scans.read_csv>`, which uses pandas to read input parameters from a csv file.
* :doc:`read_dir <../api/generated/bsmart.scans.read_dir>`, reads input files from a directory and runs tools. Useful especially for collider studies where e.g. the input files are spectrum files.
* :doc:`read_dir_mpi <../api/generated/bsmart.scans.read_dir_mpi>`, MPI version of the above, so can therefore be run over several nodes.
* :doc:`AL <../api/generated/bsmart.scans.AL>`, an active learning scan, described in the paper `Active Learning <https://arxiv.org/abs/2204.13950>`_.
* :doc:`Contour2D <../api/generated/bsmart.scans.Contour2D>`, a scan for finding points along a contour in two dimensions.
* :doc:`ContourGP <../api/generated/bsmart.scans.ContourGP>`, adapted from `excursion <https://github.com/diana-hep/excursion>`_ by Heinrich, Louppe and Cranmer, requires sklearn. Similar in aim to Contour2D, except it uses Gaussian Processes to find a contour (e.g. an exclusion curve).
* :doc:`MultiNest <../api/generated/bsmart.scans.MultiNest>`, uses `MultiNest <https://github.com/farhanferoz/MultiNest>`_. Parallelisation through MPI.
* :doc:`Diver <../api/generated/bsmart.scans.Diver>`, uses `Diver <https://diver.hepforge.org/>`_. Parallelisation through MPI.
* :doc:`MLS <../api/generated/bsmart.scans.MLS>` `Machine Learning Scan <https://arxiv.org/abs/1708.06615>`_ based on the version from `xBit <https://arxiv.org/abs/1906.03277>`_, adapted with Farid Ibrahimov.


For information about general scan settings see :mod:`bsmart.core`