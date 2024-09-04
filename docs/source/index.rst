.. ringtail documentation master file, created by
   sphinx-quickstart on Thu Apr 18 11:08:19 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
.. image:: https://user-images.githubusercontent.com/41704502/170797800-53a9d94a-932e-4936-9bea-e2d292b0c62b.png

Ringtail
##########
*A tool for handling results from virtual screening of molecules*

Ringtail is an open-source lightweight and highly customizable Python package used to organize, filter, and visualize docking data from virtual screening. Ringtail reads collections of virtual screening results in form of Docking Log File (DLG) from `AutoDock-GPU <https://github.com/ccsb-scripps/AutoDock-GPU>`_, or docking result strings or PDBQT from `AutoDock-Vina <https://github.com/ccsb-scripps/AutoDock-Vina>`_, and inserts them into an SQLite database. It then allows for the filtering of results with numerous pre-defined filtering options, generation of a simple result scatterplot, export of molecule SDFs, and export of CSVs of result data. Result file parsing is parallelized across the user's CPU (or a chosen number) to greatly enhance efficieny of the database writing.

Ringtail comes with a user-friendly :ref:`command line tool <cmdline>` for straight-forward usage, as well as an extensive :ref:`API <api>` for more advanced use in e.g., scripting. 

Ringtail is developed by the `Forli lab <https://forlilab.org/>`_ at the `Center for Computational Structural Biology (CCSB) <https://ccsb.scripps.edu>`_ at `Scripps Research <https://www.scripps.edu/>`_.

.. important:: \This publication in JCIM describes the original design, implementation, and features of Ringtail:\

      *Ringtail: A Python Tool for Efficient Management and Storage of Virtual Screening Results.*
      Althea T. Hansel-Harris, Diogo Santos-Martins, Niccolò Bruciaferri, Andreas F. Tillack, Matthew Holcomb, and Stefano Forli.
      *Journal of Chemical Information and Modeling* **2023** 63 (7), 1858-1864.
      DOI: `10.1021/acs.jcim.3c00166 <https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00166>`_

   If using Ringtail in your work, please cite this publication.


Ringtail offers a wealth of database creation and filtering options. The different sections linked below will describe each option in detail. 
To get started, first follow the instructions to :ref:`install Ringtail <installation>`, then navigate to :ref:`Getting started with Ringtail <get_started>` for a quick overview of the basic usage of Ringtail from the command line.
For more advanced and customizable use, learn how to use the :ref:`Ringtail API <api>`.

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
   database_traversing
   compare

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Full Python Documentation

   ringtail
   genindex
   modindex
