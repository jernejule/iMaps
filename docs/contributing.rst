============
Contributing
============

Installing prerequisites
========================

Make sure you have Python_ 3.6 (or higher) installed on your system.
If you don't have it yet, follow `these instructions
<https://docs.python.org/3/using/index.html>`__.

.. _Python: https://www.python.org/

Preparing environment
=====================

`Fork <https://help.github.com/articles/fork-a-repo>`__ the main
`imaps's git repository`_.

If you don't have Git installed on your system, follow `these
instructions <http://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`__.

Clone your fork (replace ``<username>`` with your GitHub account name) and
change directory::

    git clone https://github.com/<username>/imaps.git
    cd imaps

Prepare imaps for development::

    pip install -e .[docs,package,test]

.. note::

    We recommend using `pyvenv <http://docs.python.org/3/library/venv.html>`_
    to create an isolated Python environment for imaps.

.. _imaps's git repository: https://github.com/jernejule/imaps

Running tests
=============

To run all the tests, use Tox_::

    tox

To run just one test module or one specific test you can use unittest_::

    # Run module with test ``test_xyz.py``
    python -m unittest imaps.tests.test_xyz
    # Run a ``test_7`` in test case ``MyTestCase`` in  module``test_xyz.py``
    python -m unittest imaps.tests.test_xyz.MyTestCase.test_7

.. _Tox: http://tox.testrun.org/
.. _unittest:  https://docs.python.org/3/library/unittest.html

Building documentation
======================

.. code-block:: none

    python setup.py build_sphinx

Preparing release
=================

Checkout the latest code and create a release branch::

    git checkout master
    git pull
    git checkout -b release-<new-version>

Make the following changes:

    - Replace the *Unreleased* heading in ``docs/CHANGELOG.rst`` with the new
      version, followed by release's date (e.g. *13.2.0 - 2018-10-23*).
    - Bump version in the ``imaps/__about__.py`` file.

.. note::

    Use `Semantic versioning`_.

Commit changes to git::

    git commit -a -m "Prepare release <new-version>"

Push changes to your fork and open a pull request::

    git push --set-upstream <imaps-fork-name> release-<new-version>

Wait for the tests to pass and the pull request to be approved. Merge the code
to master::

    git checkout master
    git merge --ff-only release-<new-version>
    git push <imaps-upstream-name> master <new-version>

Now you can release the code on PyPI. Clean ``build`` directory::

    python setup.py clean -a

Remove previous distributions in ``dist`` directory::

    rm dist/*

Remove previous ``egg-info`` directory::

    rm -r *.egg-info

Create source distribution::

    python setup.py sdist

Build wheel::

    python setup.py bdist_wheel

Upload distribution to PyPI_::

    twine upload dist/*

.. _Semantic versioning: https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred
.. _PyPI: https://pypi.python.org/