.. ringtail documentation master file.
   Adapt freely, but keep the root toctree directives at the bottom.

.. image:: https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png
   :width: 350px

.. rst-class:: hide-title

Ringtail
##########
*A tool for handling results from virtual screening of molecules*

Ringtail is an open-source, lightweight, and highly customizable Python package for organizing, filtering, and exploring the results of molecular virtual screening — from a handful of ligands to **tens of millions**. It reads collections of docking results such as the latest SDFs from `AutoDock-6 <https://github.com/forlilab/AutoDock>`, Docking Log Files (DLGs) from `AutoDock-GPU <https://github.com/ccsb-scripps/AutoDock-GPU>`_, and PDBQTs from `AutoDock-Vina <https://github.com/ccsb-scripps/AutoDock-Vina>`_ — into a compact database that stays fast to query as it grows, backed by either DuckDB (the default) or SQLite. Result-file parsing is parallelized across your CPUs for fast database writing.

Once your docking results are in a database, Ringtail gives you a wealth of ways to apply your chemical intuition to narrow down the results to likely pharmacological hits: filter by docking score, ligand efficiency, receptor interactions, or ligand chemistry; cluster for diversity; compare hits across targets; and export exactly the molecules and data you want.

Use Ringtail however suits you — a user-friendly :ref:`command line tool <cmdline>` for straightforward screening, or an extensive :ref:`Python API <api>` for building it into your own pipelines.


Fast, even at scale
*******************

With the DuckDB backend, filtering stays in the seconds range as a library grows into the millions of ligands (DuckDB; Intel i9, 64 GB RAM, SSD):

.. csv-table::
   :header: "Ligands","Poses","Database size","Score filter","Score + interaction filter"

   "100,000","277,048","0.20 GB","1.2 s","1.4 s"
   "2,000,000","5,448,313","3.3 GB","1.4 s","5.4 s"
   "9,039,451","24,801,508","15 GB","2.0 s","16.5 s"

See :ref:`changes` for the full benchmarks, including a comparison across Ringtail versions and storage engines.


What you can do with Ringtail
*****************************

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Read any AutoDock output
      :link: get_started
      :link-type: ref

      Ingest AutoDock-GPU DLGs, AutoDock-Vina PDBQTs, and AD6 SDFs — or docking results
      streamed straight from memory, no files required.

   .. grid-item-card:: Filter the way you think
      :link: cmdline
      :link-type: ref

      Screen by docking score, ligand efficiency, percentiles, receptor interactions
      (hydrogen-bond / van der Waals / reactive), ligand substructure and 3D position, or
      molecular weight.

   .. grid-item-card:: Choose your engine
      :link: big_data
      :link-type: ref

      Store results in DuckDB (fast, the default) or SQLite — the same Ringtail, your choice
      of backend.

   .. grid-item-card:: Cluster & compare
      :link: compare
      :link-type: ref

      Reduce to diverse representatives with Morgan- or interaction-fingerprint clustering, and
      cross-reference ligands across multiple targets to find selective binders.

   .. grid-item-card:: Export what you need
      :link: api
      :link-type: ref

      Write SDFs, CSVs with exactly the columns you want (no SQL required), receptor PDBQTs, or
      a smaller standalone database.

   .. grid-item-card:: Built for big data
      :link: big_data
      :link-type: ref

      Compress databases for transfer off an HPC, merge screens run in batches, and keep working
      smoothly at tens of millions of poses.


Quick start
***********

After :ref:`installing Ringtail <installation>`, a screen is two commands:

.. code-block:: bash

   # write a folder of docking results into a database
   rt_process_vs write --docking_results results_folder/ --recursive

   # filter for the strongest binders and write a results log
   rt_process_vs read --input_db output.db --eworst -8 --output_log hits.txt

Head to :ref:`Getting started with Ringtail <get_started>` for a full walkthrough, or the
:ref:`API <api>` to script Ringtail into your own pipeline.


Citing Ringtail use
********************
Ringtail is developed by the `Forli lab <https://forlilab.org/>`_ at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_ at `Scripps Research <https://www.scripps.edu/>`_.

.. important:: \This publication in JCIM describes the original design, implementation, and features of Ringtail:\

      *Ringtail: A Python Tool for Efficient Management and Storage of Virtual Screening Results.*
      Althea T. Hansel-Harris, Diogo Santos-Martins, Niccolò Bruciaferri, Andreas F. Tillack, Matthew Holcomb, and Stefano Forli.
      *Journal of Chemical Information and Modeling* **2023** 63 (7), 1858-1864.
      DOI: `10.1021/acs.jcim.3c00166 <https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00166>`_

   If using Ringtail in your work, please cite this publication.



.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Manual

   installation
   faq
   changes
   upgrade_database

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Using Ringtail

   get_started
   cmdline
   api
   compare
   big_data
   compress

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Full Python Documentation

   modules
   genindex
   modindex
