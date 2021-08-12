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

See example <examples/usage.ipynb>`_ for example usage

Installation
------------
clone the repo.  To make a conda env with what you'll need, do the following (default env name is `mcbac-env`

.. code:: bash
   conda env create -n {optional env name} -f environment.yml

Then, to install the code do

.. code:: bash
   conda activate {env name}
   pip install -e . --no-deps



Credits
-------

This package was created with Cookiecutter_ and the `wpk-nist-gov/cookiecutter-pypackage`_ Project template forked from `audreyr/cookiecutter-pypackage`_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`wpk-nist-gov/cookiecutter-pypackage`: https://github.com/wpk-nist-gov/cookiecutter-pypackage
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
