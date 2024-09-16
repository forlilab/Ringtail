.. _installation:

Installing ringtail
###################
There are three different alternatives to installing Ringtail: through :ref:`conda-forge <condaforge>` which will install all dependencies, through the Python package manager :ref:`PyPi <pypi>` where some packages need to be installed separately, and directly from :ref:`source code <sourcecode>` for advanced users looking to make their own code changes. It is necessary to use an environment manager like conda or mamba to organize your Ringtail :ref:`environment <envsetup>` as some of the dependencies can only be installed in a managed environment. The installation instructions uses conda as an example, but you are free to use any python environment manager. Ringtail 2.0 requires Python 3.9 or higher (tested to 3.12). 

.. _pypi:
Installation from PyPI
*************************
To install Ringtail from PyPi, create then activate your :ref:`ringtail environment <envsetup>`, then simply use pip in your terminal:

.. code-block:: bash

    $ pip install ringtail

A few dependencies may be needed, including:

* meeko (another Forli lab tool)
* rdkit 
* multiprocess (only needed on MacOS)
* scipy
* pandas
* chemicalite (only available through conda-forge)

.. code-block:: bash

    $ pip install <dependency>

    $ conda install -c conda-forge chemicalite

Upgrading to a newer Ringtail version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have a previous version of Ringtail installed you can update the package by using the tag ``-U``, or by specifying the version

.. code-block:: bash

    $ pip install -U ringtail

    $ pip install ringtail==2.0

Make sure to :ref:`upgrade any databases <upgrade_database>` made with Ringtail v1 if you intend to use them with Ringtail v2.0.


.. _condaforge:
Installation from conda-forge
******************************
To install from conda-forge create a ringtail environment if needed, and run the following in the active environment:

.. code-block:: bash

    $ conda install -c conda-forge ringtail

The conda-forge installation will handle all dependencies, so no other installations are necessary. 

.. _sourcecode:
Installation from source code
******************************
To install Ringtail from source code you will need the same dependencies as for the :ref:`PyPi installation <pypi>`. 
After activating the environment, navigate to the main Ringtail ringtail directory and run:

.. code-block:: bash

    $ git clone git@github.com:forlilab/Ringtail.git
    $ cd Ringtail
    $ pip install .

This will automatically fetch the required modules and install them into the current environment.

If you wish to make the code for Ringtail **editable** without having to re-run ``pip install .``, instead use

.. code-block:: bash

    $ pip install --editable .

Test installation
------------------
If you would like to test your installation of Ringtail, or after you make changes to the code, a set of automated tests are included with the source code. To begin, you must install pytest in the Ringtail environment:

.. code-block:: bash    

    $ pip install pytest

Next, navigate to the ``test`` subdirectory within the cloned Ringtail directory and run pytest by calling

.. code-block:: bash

    $ pytest

The compounds used for the testing dataset were taken from the `NCI Diversity Set V <https://wiki.nci.nih.gov/display/NCIDTPdata/Compound+Sets>`_. The receptor used was `PDB: 4J8M <https://www.rcsb.org/structure/4J8M>`_.

.. _envsetup:
Setting up your environment
**************************
To set up your environment use for example `conda <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ or `micromamba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_, and ensure the python version is 3.9, 3.10, 3.11, or 3.12 (Ringtail 2.0.0 has not been tested for other versions). 

.. code-block:: bash

    $ conda create -n ringtail python=3.10
    $ conda activate ringtail

You can install packages from PyPi as well as other channels like ``conda-forge`` in your environment. To use PyPi/pip, you may have to first install it in your environment (especially for lightweight environment managers like micromamba). 

.. code-block:: bash

    $ conda install <package>

    $ conda install -c conda-forge <package_found_on_conda-forge>

