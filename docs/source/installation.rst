.. _installation:
Installing ringtail
==================

Installation (from PyPI)
-----------------------
Please note that Ringtail requires Python 3.9 or 3.10.

.. code-block:: bash

    $ pip install ringtail

If using conda, ``pip`` installs the package in the active environment.

Also note that if using MacOS, you may need to install Multiprocess separately:

.. code-block:: bash

    $ pip install multiprocess


Installation from source code
---------------

.. code-block:: bash

    $ conda create -n ringtail python=3.10
    $ conda activate ringtail

After activating the environment, navigate to the desired directory for installing Ringtail and do the following:

.. code-block:: bash

    $ git clone git@github.com:forlilab/Ringtail.git
    $ cd Ringtail
    $ pip install .

This will automatically fetch the required modules and install them into the current conda environment.

If you wish to make the code for Ringtail **editable** without having to re-run ``pip install .``, instead use

.. code-block:: bash

    $ pip install --editable .

Test installation
-------------------
If you would like to test your installation of Ringtail, a set of automated tests are included with the source code. To begin, you must install pytest in the Ringtail conda environment:

.. code-block:: bash    

    $ pip install -U pytest

Next, navigate to the ``test`` subdirectory within the cloned Ringtail directory and run pytest by simply calling

.. code-block:: bash

    $ pytest

The compounds used for the testing dataset were taken from the `NCI Diversity Set V <https://wiki.nci.nih.gov/display/NCIDTPdata/Compound+Sets>`_. The receptor used was `PDB: 4J8M <https://www.rcsb.org/structure/4J8M>`_.
