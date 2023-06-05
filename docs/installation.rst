Installation
===============

.. todo::
    Update these commands for installation

It is recommended to set up a separate virtual or conda environment to install and immunopepper.

The basic ImmunoPepper package can be installed using pip:

.. code-block:: console

    pip install immunopepper

Alternatively, ImmunoPepper can also be installed from source using:

.. code-block:: console

    pip install -r requirements.txt -r requirements_dev.txt
    make install

After installation, please consult the help screen for further usage options:

.. code-block:: console

    immunopepper -h

Prerequisites
-------------

ImmunoPepper takes a splicing graph as input. This splicing graph has to be generated using the
SplAdder pipeline. Further information about SplAdder is available on its `GitHub
page <https://github.com/ratschlab/spladder>`_ or the `Online
documentation <https://spladder.readthedocs.io/en/latest/>`_.






