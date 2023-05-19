#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file.

import os
import sys

from setuptools import setup

# To check version before setup, run:
# python setup.py --version

# To push to PyPi:
# python setup.py sdist bdist_wheel --universal
# twine upload dist/*  
# python -m twine upload dist/*  

# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    tox -e test

If you don't already have tox installed, you can install it with:

    pip install tox

If you only want to run part of the test suite, you can also use pytest
directly with::

    pip install -e .[test]
    pytest

For more information, see:

  http://docs.astropy.org/en/latest/development/testguide.html#running-tests
"""

if 'test' in sys.argv:
    print(TEST_HELP)
    sys.exit(1)

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py build_docs'. Instead you will need to run:

    tox -e build_docs

If you don't already have tox installed, you can install it with:

    pip install tox

You can also build the documentation with Sphinx directly using::

    pip install -e .[docs]
    cd docs
    make html

For more information, see:

  http://docs.astropy.org/en/latest/install.html#builddocs
"""

if 'build_docs' in sys.argv or 'build_sphinx' in sys.argv:
    print(DOCS_HELP)
    sys.exit(1)

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    version = get_version(root='..', relative_to=__file__)
except Exception:
    version = '{version}'
""".lstrip()

setup(
    version='2.0.0.dev2'
    # use_scm_version={
    #     'write_to': os.path.join('dorado', 'version.py'),
    #     'write_to_template': VERSION_TEMPLATE,
    #     "local_scheme": "no-local-version",
    #     "version_scheme": "python-simplified-semver"
    #                    }
        )
# where is this file? os.path.dirname(__file__)
# import  astropy.config as acfg
# import ./dorado/config


# Entry points and automatic script creation
# Setuptools supports automatic creation of scripts upon installation, that run code within your package if you specify them as entry points. An example of how this feature can be used in pip: it allows you to run commands like pip install instead of having to type python -m pip install.

# https://setuptools.pypa.io/en/latest/userguide/quickstart.html






