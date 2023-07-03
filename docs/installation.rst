Installation
===============

Installation from source
------------------------

.. note::
    For the installation from source, the user needs to clone the git repository of immunopepper.
    This can be done by running the following command:

    .. code-block:: console

        git clone https://github.com/ratschlab/immunopepper.git


It is recommended to set up a separate virtual or conda environment to install immunopepper.

**The python version in the environment should be >= 3.9**

.. code-block::

        conda create -n immunopepper python=3.9
        conda activate immunopepper

The user should also check if *gcc* is installed on the system, as it is required for the installation of the tool.
This can be checked by running the following command:

.. code-block::

    gcc --version


The rest of the installation can be performed by running:

.. code-block::

    conda install cython
    conda install -c bioconda 'pyvcf3==1.0.3'
    make install


After installation, please consult the help screen for further usage options:

.. code-block::

    immunopepper -h


Prerequisites
-------------

ImmunoPepper takes a splicing graph as input. This splicing graph has to be generated using the
SplAdder pipeline. Further information about SplAdder is available on its `GitHub
page <https://github.com/ratschlab/spladder>`_ or the `Online
documentation <https://spladder.readthedocs.io/en/latest/>`_.

Mode cancerspecif
-----------------

In this mode Java needs to be installed. The user can check if Java is installed by running the following command:

.. code-block::

    java -version

Mode mhcbind
------------

If the mode :ref:`mhcbind <mhcbind>` wants to be run  , the user needs to clone the mhcbind github repository inside the main immunopepper folder.

This can be done by running the following command:

.. code-block::

    git clone https://github.com/openvax/mhctools.git

.. todo::
    Check if this is done like this



