=====
mcbac
=====


.. image:: https://img.shields.io/pypi/v/mcbac.svg
        :target: https://pypi.python.org/pypi/mcbac

.. image:: https://img.shields.io/travis/wpk-nist-gov/mcbac.svg
        :target: https://travis-ci.com/wpk-nist-gov/mcbac

.. image:: https://readthedocs.org/projects/mcbac/badge/?version=latest
        :target: https://mcbac.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




mCBAC implementation


* Free software: NIST license
* Documentation: https://mcbac.readthedocs.io.


Features
--------

Python version of m-CBAC calculation (see JiLink_)

See `example <examples/usage.ipynb>`_ for example usage

Installation
------------
clone the repo.  To make a conda env with what you'll need, do the following (default env name is `mcbac-env`

.. code-block:: console

   $ conda env create -n {optional env name} -f environment.yml

To optionally install development packages:

.. code-block:: console

   $ conda env update -n {env name} -f environment-dev.yml

Then, to install the code do

.. code-block:: console

   $ conda activate {env name}
   $ pip install -e . --no-deps


To run tests, use

.. code-block:: console

   # tests on local env
   $ pytest

   # tests on all python versions
   $ tox


Optionally, to install directly from the repo via pip, do the following:

.. code-block:: console

   $ pip install git+https://github.com/wpk-nist-gov/mcbac.git@develop







Credits
-------

This package was created with Cookiecutter_ and the `wpk-nist-gov/cookiecutter-pypackage`_ Project template forked from `audreyr/cookiecutter-pypackage`_.


.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`wpk-nist-gov/cookiecutter-pypackage`: https://github.com/wpk-nist-gov/cookiecutter-pypackage
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _JiLink: https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.0c01524
